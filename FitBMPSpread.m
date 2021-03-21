function FitBMPSpread
% applies to localized BMP source and its distribution in space either as a fn
% of time or of applied Noggin.
%
% Notes on data sets:
%   '*Basal*_250ng*' the fixed threshold for detecting BMP cells is a bit
% too high for a few nuclei. In for several BMP secreting cells with apical noggin, one
% finds low pSmad1 in a few nominally secreting cells and their neighbors. This
% is evident in the histogram of pSmad1 within the BMP mask. Presumably the
% DOX regulation is off, though the nuclear marker is still strong.
%   '*Apical*' a semilog plot in hist_all_keys is quite close to linear.
data_dir = '/Users/siggia/Dropbox/Desktop/phan/image_proc/BMPpropagation/202003_BMP4Propagation_4h_8h_12h_24h/';
all_keys = {'*8h*', '*12h*', '*24h*'};

data_dir = '/Users/siggia/Desktop/phan/image_proc/Noggin_apical_basal/202003_NOGGINsApicalBasalSensitivity/';
all_keys = {'*Apical*_10ng*', '*Apical*_50ng*', '*Apical*_250ng*'};
% all_keys = {'*Basal*_10ng*', '*Basal*_50ng*', '*Basal*_250ng*'};

% merging here inlet, mid, outlet. For the low flow case, data variable as
% regards smad1 cloud, high flow more uniform. Doubled bin size in
% hist_all_keys to get better plots, Decay lengths of fits 19,20mu lo/hi but levels
% a bit different
data_dir = '/Users/siggia/Desktop/phan/image_proc/BMP500_1flowcell/20200218_WTandBMP_500to1_20200213-selected/';
all_keys = {'*Low*', '*High*'};

% %case of noggin spreading and inhibiting uniform BMP.  Does not work easily
% %since need invert ON/OFF for smad1 and issues around boundaries where
% %noggin cells not imaged but still see inhibition. Also problems with STD
% %since only 1 image given
% data_dir = '/Volumes/Data/tien_polarity_paper/202008_NOGGIN_SPREAD/';
% all_keys = {'*24h*'};

% averages over all the images corresponding to one key 
for i = 1:length(all_keys)
    stats(i) = fit1data_type(data_dir, all_keys{i});
end

% plots pSmad vs distance for all the keys in one graph. 
hist_all_keys(stats, all_keys);
return

function stats_out = fit1data_type(data_dir, key)
%
% pull out all files in data_dir with reg-ex pattern=key (eg "*24h*"), and
% return composite statistics with 

verbose = 1;
r_nuc = 8;   % determined from ImageJ by counting pixels
med_filt = 5;
dst_step = 8;    % collect stats vs distance in these bins.
dst_mx = 192;   % max dst from source to collect pSmad1, < size_img/2
frac_bmp = 1/200;
idapi=1;  ibmp=2;  ismad=3;
% channels DAPI, BMP, pSMAD, and 12 bit images
all_img = dir([data_dir, key]);
nimg = length(all_img);
if nimg < 2
    fprintf('< 2 files found matching key= %s, bye!\n', key);
    return
end

for i = 1:nimg
    img = read_img([data_dir,all_img(i).name]);
    nchannels = size(img, 3);
    for j = 1:nchannels
        img(:,:,j) = medfilt2(img(:,:,j),[med_filt, med_filt]);
    end
    fprintf('\nread file= %s, max= %d %d %d after median flt\n',all_img(i).name, squeeze(max(max(img))) );
    img = uint8((img+1)/16);
    bmp_mask = img(:,:,ibmp);
    % based on fraction of bmp, can still be badly off ie too low if only a
    % few bmp cells in image by chance.
    thresh_q = quantile(bmp_mask(:), (1-frac_bmp)); 
    % for BMP500_1flowcell data, adjust literal thresh, different microscope
    % than other data.
    if strcmp(key, '*Low*')
        thresh = max(thresh_q, 8);
        thresh = min(thresh, 35);
    elseif strcmp(key, '*High*')
        thresh = max(thresh_q, 14);
    else
        thresh = 45;
    end
    fprintf('threshold quantile= %d, using %d\n', thresh_q, thresh);
    bmp_mask = bmp_mask > thresh;
    bmp_dst = bwdist(bmp_mask);
    if verbose
        check_mask(img(:,:,ibmp), bmp_mask);
        title(sprintf('BMP mask for thresh= %d, fname= %s',thresh, key ));
  
        % histogram of smad1 within BMP mask:
        data = img(:,:,ismad);
        data = data(bmp_dst == 0);
        figure, hist(double(data), 32);
        title(sprintf('pSmad within BMP mask. pts= %d, fname= %s',numel(data), key) );
    end
    % background value of pSmad1 in nuclei far from source,assuming >~ 1/2
    % area is nuclei and fraction of nuclei with pSmad1 ON is < 1/2
    dapi_mask = img(:,:,idapi);
    thresh_dapi = quantile(dapi_mask(:), 0.75);
    dapi_mask = img(:,:,idapi) >= thresh_dapi;
    data = img(:,:,ismad);
    data = data(dapi_mask);
    thresh_smad = quantile(data, 0.10);
    fprintf('background level pSmad1 in nuclei= %d, (subtr from pSmad) based on dapi_thresh= %d\n', thresh_smad, thresh_dapi);
    % collect stats by dst from BMP 
    stats(i) = count_by_dst(img(:,:,ismad), bmp_dst, dst_step, dst_mx);
    stats(i).backgnd = thresh_smad;
    stats(i).fname = all_img(i).name;
end
    
stats_out = hist_stats(stats);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mask = inside(img, dst2b)
% return a mask with all points >= dst2b from boundaries
[m,n]=size(img);
mask = true(m,n);
mask(1:dst2b,:) = false;
mask((m-dst2b+1):m,:) = false;
mask(:,1:dst2b) = false;
mask(:,(n-dst2b+1):n) = false;
return

function stats = count_by_dst(data, dst_map, step, dst_mx)
% given distance matrix from source, histogram data as fn of distance,
% omitting points around edge.

nbins = floor((dst_mx+1)/step);  % assuming lims(1) = 0
lims = (0:nbins)*step;    % lims(i) <= bin(i) < lims(i+1)
npts = zeros(1,nbins);
cnts = zeros(1,nbins);
for i = 1:nbins
    mask = (dst_map >= lims(i))&(dst_map < lims(i+1))&inside(data, i*step-1);
    npts(i) = nnz(mask);
    cnts(i) = sum(data(mask));
end
[~, idx] = max(npts);
vals = (idx-1):(idx+1);
vals = max(vals, 1);  vals = min(vals, nbins);
% when area of mask starts to decrease, hint to stop binning.
fprintf('count_by_dst: max pts bin= [%d, %d) vals= %d %d %d \n',...
    lims(idx),lims(idx+1), npts(vals) );
stats.lims = lims;  stats.npts = npts; stats.cnts = cnts; 
stats.backgnd = [];  stats.fname = [];  % place holders
return

function stats_out = hist_stats(stats)
% average the stats from all images of a given type. 
% return composite stats, for plotting with data of different types.
nimg = length(stats);
lims = stats(1).lims;
nbins = length(lims)-1;

% extract npts(nimg, nbins) arrays.
npts = reshape([stats.npts], [], nimg)' + 1;  % + 1 for zero npts
cnts = reshape([stats.cnts], [], nimg)';

% remove background
for i = 1:nimg
    cnts(i,:) = max(cnts(i,:) - double(stats(i).backgnd)*npts(i,:), 0) ;
end
mean_all = cnts ./ npts;

step = floor(nbins/10);
fprintf('sample of means by img for subset of bins j= 1,%d,%d..: before removing backgnd=...\n', step, 2*step)
disp([stats.backgnd]);
for i = 1:nimg
    for j = 1:step:min(nbins, 10*step+1)
        fprintf('%5.1f  ', mean_all(i,j));
    end
    fprintf('\n');
end

stats_out.cnts = sum(cnts,1);
stats_out.std  = std(cnts, 0, 1);    % std over nimg, normed by N-1
stats_out.npts = sum(npts,1);
stats_out.lims = stats(1).lims;
figure, plot(lims(2:end), stats_out.npts)
title('total number of points in bins defined by dst from BMP');

mean = (stats_out.cnts) ./ (stats_out.npts);
bins = floor( (lims(1:nbins) + lims(2:(nbins+1)))/2 );
figure, bar(bins, mean);
title(sprintf('mean pSmad1 vs dst to BMP, for files.. %s',stats(1).fname) );
return


