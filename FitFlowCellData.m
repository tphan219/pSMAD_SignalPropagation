function scaled_smad = FitFlowCellData(verbose)

% collect the directories, file names and file contents for various conditions
% run this function multiple times to collect scaled_smad distributions and then compare
% via KS test, in compare_smad_dst().
% NB besure to set parameters img() and frac_cells().
%
% verbose = 1       do graphics and diagnostics from stats_smad_dst()
%
% see also  stats_smad_dst,  compare_smad_dst,  flowcell_dapi2colony_mask
% there are buried parameters in stats*dst and flow*mask
%

% data_dir = '/Users/siggia/Desktop/phan/image_proc/flow_cell_data/';
% fn10inlet = '20191115_MAX_Stack_exp20191115_WT BMP4 NOGGIN_100 to 20  to 1_flowRate10ulph_Inlet copy.tif';
% fn10outlet = '20191115_MAX_Stack_exp20191115_WT BMP4 NOGGIN_100 to 20  to 1_flowRate10ulph_Outlet copy.tif';
% img = struct('type', {'dapi', 'bmp', 'noggin', 'psmad'});
% frac_cells = [1, 0.05, 0.01, 1];  % following 'type' the number fraction of each stain
% fn2read = fn10outlet;

% data_dir = '/Volumes/Data 1/';   % at home
% data_dir = ['/Volumes/archive/202002_Tien_MicrofluidicData_MaximumProjection/',...
% '20200220_MAX_Stack_exp20191114_ BMP4 NOGGIN_300  to 1 _DAPI_NOGGINcells_pSMAD/'];
% fn10inlet = '20200220_MAX_Stack_exp20191114_ BMP4 NOGGIN_300  to 1_flowRate10ulph_Inlet.tif';
% fn10outlet = '20200220_MAX_Stack_exp20191114_ BMP4 NOGGIN_300  to 1_flowRate10ulph_Outlet.tif';
% fn1000inlet = '20200220_MAX_Stack_exp20191114_ BMP4 NOGGIN_300  to 1_flowRate1000ulph_Inlet.tif';
% fn1000outlet = '20200220_MAX_Stack_exp20191114_ BMP4 NOGGIN_300  to 1_flowRate1000ulph_Outlet.tif';
% img = struct('type', {'dapi', 'noggin', 'psmad'});
% frac_cells = [1, 1/300, 1];
% % results: compared for this DS, outlet to inlet  kstest2:  prob= 9.77e-09, kstest= 2.92e-01

% % the third tifx(:,:,1) from this series == sum of dapi + smad, and (:,:,2:3) == dapi,
% % ignore. 
% data_dir = ['/Volumes/archive/202002_Tien_MicrofluidicData_MaximumProjection/',...
% '20191025_WT and BMP4_10 to 1_exp20190927_DAPI_pSMAD_Merged/'];
% fn10  = '20191025_WT and BMP4_10 to 1_flowRate10ulph_exp20190927_MAX_Stack.tif';
% fn800 = '20191025_WTandBMP4_10 to 1_flowRate800ulph_exp20190927_MAX_Stack.tif';
% img = struct('type', {'dapi', 'psmad'});

% This is data used in Fig 4 of paper, and plot of Smad vs dst to noggin
% for high flow, placed in flow_cell_data/
data_dir='/Volumes/Data/tien_polarity_paper/20200311_WTandBMPandNOGGIN_300to30to1_withReporter/';
fn1000outlet = 'MAX_HighFlow1mlph_Outlet_4MpmL_20hculture_6hDox_WTBMP4NOGGIN_300to30to1_20200313_2020_03_13__18_51_57.tif';
fn1000inlet = 'MAX_HighFlow1mlph_Inlet_4MpmL_20hculture_6hDox_WTBMP4NOGGIN_300to30to1_20200313_2020_03_13__11_50_05.tif';
% data_dir = ['/Volumes/archive/202002_Tien_MicrofluidicData_MaximumProjection/',...
%     '20200311_WTandBMP_10to1_withReporter/'];
% fn1000inlet = 'MAX_HighFlow 1mLph_Inlet_4MpmL_32hculture_6hDox_WTBMP4_10to1_20200311_DAPI_H2BmCitrine_pSMAD15.tif';
% fn1000outlet = 'MAX_HighFlow 1mLph_Outlet_4MpmL_32hculture_6hDox_WTBMP4_10to1_20200311_DAPI_H2BmCitrine_pSMAD15.tif';
% fn10inlet = 'MAX_LowFlow10ulph_Intlet_4MpmL_32hculture_6hDox_WTBMP4_10to1_20200311_DAPI_H2BmCitrine_pSMAD15.tif';
% fn10outlet = 'MAX_LowFlow10ulph_Outlet_4MpmL_32hculture_6hDox_WTBMP4_10to1_20200311_DAPI_H2BmCitrine_pSMAD15.tif';
img = struct('type', {'dapi', 'bmp', 'noggin', 'psmad'});
frac_cells = [1, 1/10,1/300, 1];
fn2read = fn1000outlet;
% fnout = '20200220flow1000inlet.tif';  % not currently used

% For the localized Noggin, uniform basal BMP on filter various times. Last data run.
% NOT FLOW CELL DATA.  Need to fiddle flowcell*mask since PDF of DAPI had
% no dip to distinguish in/out colony mask. stats*dst has parameters for
% spacing of histogram etc.
data_dir = '/Volumes/Data/tien_polarity_paper/202008_NOGGIN_SPREAD/';
% fn2read = 'MAX_4h.tif';
fn2read = 'MAX_8h.tif';
fn2read = 'MAX_8h_take2.tif';
% fn2read = 'MAX_14h.tif';
% fn2read = 'MAX_24h.tif';
% fn2read = 'MAX_24h_take2.tif';
img = struct('type', {'dapi', 'noggin', 'psmad'});
frac_cells = [1, 1/300, 1];

% NB test for exact match to {'dapi', 'bmp', 'noggin', 'psmad'}, so use these names

% exclude case where raw data has overlay of primary channels as last layer.
nimg = numel(imfinfo([data_dir, fn2read]) );
if nimg ~= length(img)
    fprintf('FitFlowCellData: WARNING num-img from fn= %d, num-img-fields= %d, taking min\n', nimg, length(img));
    nimg = min(nimg, length(img));
end
for i = 1:nimg
    tifx = imread([data_dir, fn2read], i);
    img(i).data = tifx(:,:,1);  % seems to read RGB for what is gray scale image. For some data.
    clear tifx;
    if isa(img(i).data, 'uint16') 
        lo_hi = stretchlim(img(i).data, [0.0, 0.99]);
        img(i).data = uint8(imadjust(img(i).data, lo_hi, [0,1])/256 );  % should (..  +1)/256
        fprintf('rescaled input data to uint8 from unint16 to use [0,255]\n')
    elseif isa(img(i).data, 'uint8')
    else
        disp('FitFlowCellData: unrecognized data type, bye!');
    end
        
end
fprintf('FitFlowCellData: read %d images, each of size= %d %d from %s\n',nimg, size(img(1).data), fn2read);

idapi = find_name({img.type}, 'dapi');
ibmp = find_name({img.type}, 'bmp');
inog = find_name({img.type}, 'noggin');
ismad = find_name({img.type}, 'psmad');
fields2binary = sort([ibmp, inog]);

% With noggin data, eliminate saturating smad pixels, which are generally garbage. When
% smad only output global histogram and can trim saturated values then 
if ~isempty(inog)  % isempty(strfind(fnread, '20191025_WT and BMP4_10 to 1_exp20190927_DAPI_pSMAD_Merged') )
    disp('FitFlowCellData: resetting saturated smad pixels to mean');
    img(ismad).data(img(ismad).data==255) = mean(img(ismad).data(:));
end

colony = flowcell_dapi2colony_mask(fn2read, img(idapi).data);

% RGB projected image for ImageJ (could put colony mask as white.)
if isempty(fields2binary)  % compute outline of colony in last channel??
    imwrite(cat(3, img(idapi).data, img(ismad).data, uint8(zeros(size(img(ismad).data))) ), 'all_data.tif');
elseif ~isempty(inog)
    imwrite(cat(3, img(idapi).data, img(inog).data, img(ismad).data), 'all_data.tif');
else
    fprintf('FitFlowCellData: did not write ImageJ color image, %d %d %d %d\n',...
        idapi, ibmp, inog, ismad);
end

if verbose
    check_mask(img(ismad).data, colony)
    title(['test of colony wide mask based on DAPI image ', fn2read(1:min(12,end))]);
end

% threshold the fields corresponding to bmp and noggin secreting cells and
% make masks showing the cells that are on.
for i = fields2binary
   img(i).data(~colony) = 0;  % or replace by median so as not to disrupt thresholding??
   [cutoff, img(i).mask, areas] = thresh_data(img(i).data, frac_cells(i));
   fprintf('FitFlowCellData: thresholded= %d for %s data, #clumps= %d, avr-area= %d, frac-colony= %7.5f\n',...
       cutoff, img(i).type, length(areas), mean(areas), sum(areas)/sum(colony(:)) );
%     npts0 = sum(img(i).mask(:));
%     img(i).mask = imopen(img(i).mask, strel('disk',1));
%     npts1 = sum(img(i).mask(:));
%     fprintf('binarize data= %s, otsu thres= %d, pct mask= %5.3f, post imopen= %5.3f\n',...
%         img(i).type, cutoff, 100*[npts0, npts1]/numel(img(i).mask(:)) );
    if verbose
        check_mask(img(i).data, img(i).mask);
        title(['test image mask for data= ', img(i).type, ' ',fn2read(1:min(12,end))]);
    end
end

% when noggin data, overlay noggin islands on smad, and return stats of smad intensity
% relative to noggin islands.  Can also process distance from BMP
if ~isempty(inog)
    if verbose
        check_mask(img(ismad).data, img(inog).mask);
        title(['noggin mask superimposed on smad data ', fn2read(1:min(12,end))])
    end
%     scaled_smad = stats_smad_dst(img(ibmp), img(ismad), colony, verbose);
    scaled_smad = stats_smad_dst(img(inog), img(ismad), colony, fn2read, verbose);
else % smad1 only no noggin
    data = img(ismad).data(colony);
    [cnts, edges] = histcounts(double(data), (0:256)-0.5);
    figure, bar(0:255, cnts);
    scaled_smad = cnts;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%


function [cutoff, mask, areas] = thresh_data(img, frac_on)
% use expected fraction of active cells to guide thresholding. Return areas of spots that
% pass thresholds.

min_area = 20;
% can gaussian filter radius 1, but lowers computed thresh, need incr area to get rid of
% small weak blobs. 
cutoff1 = round(255*graythresh(img) );
cutoff2 = min(quantile(img(:), 1-frac_on), 254);  % Using mask > cutoff so if 255 would get []
cutoff = max(cutoff1, cutoff2);

mask = img > cutoff;
mask = imopen(mask, strel('disk', max(1, round(sqrt(min_area/pi))) ));  % get rid of stringy things
cc = bwconncomp(mask);
areas = regionprops(cc, 'Area');
areas = [areas.Area];
intensity = zeros(size(areas));
for i = 1:length(areas)   % get rid of remaining isolated bits  %% imopen takes care
    intensity(i) = mean(img(cc.PixelIdxList{i}));
    if length(cc.PixelIdxList{i}) >= min_area
        continue
    end
    mask(cc.PixelIdxList{i}) = 0;
end
thresh_area = quantile(areas(areas>=min_area), 0.5);
thresh_inten = max(mean(intensity), 2*cutoff);
figure, plot(min(areas, thresh_area), min(intensity, double(thresh_inten)), '.');
title(sprintf('mean intensity vs area of all cc of mask, pre filtering. Area, mean cutoffs= %d, %d',...
    thresh_area, thresh_inten));
fprintf('For mask: cutoff= max(%4.1f, %4.1f), frag (< %d) discarded= %d, retained= %d\n',...
cutoff1, cutoff2, min_area, sum(areas < min_area), sum(areas >= min_area) );
areas = areas(areas >= min_area);
return

function ii = find_name(cell_array, pattern)
% find location of exact match to pattern in cell_array, or print warning if not found
ii = find(strcmp(cell_array, pattern));
if isempty(ii)
    fprintf('could not match %s in cell array\n', pattern);
end