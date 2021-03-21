function colony = flowcell_dapi2colony_mask(fn2read, img)
% convert uint4 image to binary mask, with specific instructions for certain data sets as
% defined in fn2read

% find the outline of the DAPI image and locate any significant holes, and
% remove rough edges and odd appendages. 
if ~isa(img, 'uint8')
    fprintf('flowcell_dapi2colony_mask: img is not unit8, bye')
    colony = [];
    return
end
nuc_rad = 10;
min_area_frac = 0.99;  % insist at least this fraction of area be inside colony
colony = imgaussfilt(img, 2*nuc_rad);  % nuclear diameter around 8
cnts = imhist(colony);
cnts = filter_hist(cnts);
% find local min, which is min between backgroud peak and real DAPI.
dd_cnts = cnts(3:256) + cnts(1:254) - 2*cnts(2:255);
loc_min = (cnts(3:256)-cnts(2:255)).*(cnts(2:255)-cnts(1:254));
idx = find( (dd_cnts > 0) & (loc_min < 0));
[~, imx] = max(cnts(idx));
imx = idx(imx);   % index into dd_cnts array with max value of cnts if ambiguous.
imx = imx + 1;
area_frac = cumsum(cnts)/sum(cnts);
imx2 = find(area_frac> 1-min_area_frac,1);
fprintf('flowcell_dapi2colony_mask: colony >= Min(%d based on dip in PDF, =%d for min_area= %4.3f)\n',...
    imx, imx2, min_area_frac);
imx = min(imx, imx2);

figure, bar(0:255, cnts)
title(sprintf('hist gaussfilt dapi, radius= %d, mask== (dapi >= %d), pct area out= %3.1f',...
    2*nuc_rad, imx, area_frac(imx)));
colony = colony >= imx;

% colony = imopen(colony, strel('disk', 2*nuc_rad)); % there are isolated cells to eliminate, handled by largest_cc()
colony = imclose(colony, strel('disk', 2*nuc_rad));  % small holes to fill
if strcmp(fn2read, '20191115_MAX_Stack_exp20191115_WT BMP4 NOGGIN_100 to 20  to 1_flowRate10ulph_Outlet') 
    % case filamentacious stuff
    radius = 80;
    fprintf('doing erode, then dilate radius= %d around conn comp, fn= %s\n',radius, fn2read); 
    str_disk = strel('disk', radius);
    colony = imerode(colony, str_disk);
    colony = largest_cc(colony);
    colony = imdilate(colony, str_disk);
elseif strcmp(fn2read, '20200220_MAX_Stack_exp20191114_ BMP4 NOGGIN_300  to 1_flowRate1000ulph_Inlet.tif')
     % 20 columns of 255, 0 at end of DAPI img, infects all images.
    colony(:, (end-19):end) = 0;
else
    colony = largest_cc(colony);
end
colony = imerode(colony, strel('disk', 2*nuc_rad));
fprintf('flowcell_smad2colony_mask: largest cc= %d, area image= %d (%d x %d), final colony= %d\n',...
    sum(colony(:)), numel(colony), size(colony), sum(colony(:)) );

%%%%%%%%%%%%%%%%%%
function colony = largest_cc(colony)
% return largest cc and 
cc = bwconncomp(colony);
areas = regionprops(cc, 'Area');
areas = [areas.Area];
[~, imx] = max(areas);
for i = 1:length(areas)   % get rid of remaining isolated bits
    if i==imx
        continue
    end
    colony(cc.PixelIdxList{i}) = 0;
end

function cnts = filter_hist(cnts)
% smooth histogram of uint8 data cnts(1:256)
sig = 10;
cnts = reshape(cnts, 1, 256);
cnts = fft(cnts);
ksq = 1:256;
ksq = min([ (ksq-1).^2; (257-ksq).^2]);
fac = exp(-ksq/(2*sig^2));
cnts = cnts .* fac;
cnts = ifft(cnts);
return