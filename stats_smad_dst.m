function scaled_smad = stats_smad_dst(marker, smad, colony, fname, verbose)
%
% For FitFlowCellData: overlay noggin and pSmad images, bin pSmad in rings defined by
% bwdist function, but contained in square centered around each noggin source that is just
% large enough to include entire ring. NB distance rings with respect to all Noggin source
% as well as exterior of colony.  
%`  Returned statistics based on quantile's ie level of smad such that ~10% of pixel are
% over this level. Avoids problems with bad cells that can saturate the smad1 signal at
% 255, and quantile ~90% restricts to mostly nuclei. To eliminate large scale variation in images, 
% report max_r(smad(r)) - smad(r = 1)
%
% input:   (NB box_side2 parameter to reset, BMP vs NOGGIN in marker)
%   marker      struct field for type, data, mask
%   smad        struct with field type data
%   colony      binary, define colony
%   verbose     print stats for each noggin center
%
% output:
%   scaled_smad     (n-radii-bin, noggin-centers) array quantile thresholds

% average smad signal as function of distance from nearest Noggin cell.
% compute average in uniformly spaced distance distance-radius bins 
nbin = 20;
box_side2 = 15;  % spacing of radii bins use ~3 for BMP source turning on Smad, ~20 for noggin inhibition

% treat outside of colony as noggin on, thus dst rings terminate that dst from colony edge
dst = bwdist(marker.mask | ~colony); 
% dst = bwdist(marker.mask);
% dst(~colony) = 0;   

% over = uint8(dst);  %% irrelevant image
% figure, imshow(cat(3,over+ 255*uint8(marker.mask), over, over));
% title('dst from noggin cells (red)')
 
cc = bwconncomp(marker.mask);
props = regionprops(cc, 'basic');
areas = [props.Area];
[areas, isort] = sort(areas, 'descend');
props = props(isort);

% find max side of bounding box around noggin center.
side_nog = reshape([props.BoundingBox], 4, []);
side_nog(1:2,:) = [];   % remove centroid
side_nog = ceil(max(side_nog)/2 );
% mean smad, nog calculated for the noggin ON mask centers
mean_nog = zeros(1, cc.NumObjects);
mean_smad = zeros(1, cc.NumObjects);
for i = 1:cc.NumObjects 
    mean_nog(i) = mean(marker.data(cc.PixelIdxList{i}));
    mean_smad(i) = mean(smad.data(cc.PixelIdxList{i}));
end

quantile_cutoff = 0.9;
all_smad = smad.data(colony);
fprintf('stats_smad_dst: mean pSmad(colony)= %4.1f; for quantile= %4.2f, cutoff= %3.1f, #noggin= %d\n',...
    mean(double(all_smad)), quantile_cutoff, quantile(double(all_smad), quantile_cutoff), cc.NumObjects );
scaled_smad = zeros(nbin, cc.NumObjects);
for i = (1:nbin)
    shell = box_side2*(i-1) < dst & dst <= box_side2*i;
    toprint = 200;  % if verbose then number of noggin centers to print.
    mod2print = max(floor(cc.NumObjects/toprint), 1);
    for j = 1:cc.NumObjects
        ctr = round(props(j).Centroid);
        % cutout a square large enough to include everything in shell. NB is dst from nearest
        % mask pixel, while box computed from centroid. NB ctr input as xy, cut_out flips
        % to ij. Create a mask based on cut out on shell,then extract data within mask.
        pixels = cut_out(ctr, i*box_side2+side_nog(j), shell);
        if isempty(pixels)
            thresh = 0;
        else
            data = cut_out(ctr, i*box_side2+side_nog(j), smad.data);
            data = double(data(pixels));
            thresh = quantile(data, quantile_cutoff); % set quantile, 90% better contrast than 60%
        end
        scaled_smad(i, j) = thresh;  % mean(data(data>=thresh) ); % thresh avoids a few cases with 1 nucl==255, now set to median
        if verbose && i==1 && mod(j,mod2print)==0
            fprintf('stats_smad_dst: noggin mask: j= %d, ctr= %d %d, area= %d, mean(nog,smad)= %3.0f, %3.0f, npts>quantile= %4.1f \n',...
                j, ctr, props(j).Area, mean_nog(j), mean_smad(j), sum(pixels(:))*(1-quantile_cutoff) );
        end
    end
end

% stats
% mostly contrast of quartile vs mean, ring=1 and on noggin secreting cell
% should both be smad1 OFF.
figure, plot(mean_smad, scaled_smad(1,:), '*')
title('quartile smad ring 1 vs mean smad in noggin region');

% take max over all rings/bins and compare with first ring that --> 0
[mx, ix] = max(scaled_smad);
figure, histogram(mx - scaled_smad(1,:), 15);
title(sprintf('max smad1 ring - ctr. nbin= %d, dlt-R= %d, pct keep= %d, #nog= %d',...
    nbin, box_side2, round(100*(1-quantile_cutoff)), cc.NumObjects ));

% scatter plot for deviation from max vs distance/ring for each center.
figure, plot(ix, mx - scaled_smad(1,:), '*');
title(sprintf('max dlt-smad1 vs bin. nbin= %d, dlt-R= %d, pct keep= %d, #nog= %d',...
    nbin, box_side2, round(100*(1-quantile_cutoff)), cc.NumObjects ));

% plot average of Smad and std of mean vs dst from each noggin source.
smad2plot = scaled_smad - ones(nbin,1)*scaled_smad(1,:);
smad2plot = max(smad2plot, 0);
scl = 425.0961/1024;  %% for flowcell data, microns/pixel, same Noggin inhibition
dst = scl*(1:nbin)*box_side2; % scl*((1:nbin)+0.5)*box_side2;
figure
% plot(dst, mean(smad2plot,2), '*k');
errorbar(dst, mean(smad2plot,2), std(smad2plot,0,2)/sqrt(size(scaled_smad,2)), '*k');
xlabel('Distance to closest Noggin secreting cell ({\mu}m)', 'FontSize',18);
title([fname(1:min(12,end)),' pSMAD1 quartile of ring relative to ctr with STD of average'])

mtx2save = zeros(3,length(dst));
mtx2save(1,:) = dst;
mtx2save(2,:) = mean(smad2plot,2);
mtx2save(3,:) = std(smad2plot,0,2)/sqrt(size(scaled_smad,2));
csvwrite([fname, '.csv'], mtx2save');
return

%%%%%%%%%%%%%%%%%
function pixels = pixels_in_shell(size_mtx, radii)
% for list of radii find all pixel, in shells.  NOT CHECKED
for n = 1:(length(radii)-1)
    ijpix = zeros(2,2*radii(n+1));
    ptr = 0;
    rng = (-radii(n+1)):radii(n+1);
    for i = rng
        for j = rng
            dst = i^2 + j^2;
            if( radii(n)^2 < dst && dst <= radii(n+1)^2)
                ptr = ptr + 1;
                ijpix(:, ptr) = [i; j];
            end
        end
    end
    ijpix(:, (ptr+1):end) = [];
    ijpix = ijpix + size_mtx' * ones(1, size(ijpix, 2));
    pixels{n} = sub2idn(size_mtx, ijpix(1,:), ijpix(2,:));
end
        
function cutout = cut_out(ctr, side2, data)
[m,n] = size(data);
cutout =  data(max(1, ctr(2)-side2):min(ctr(2)+side2, m), max(1, ctr(1)-side2):min(ctr(1)+side2, n) );
