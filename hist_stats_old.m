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

