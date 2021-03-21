function hist_all_keys(stats, keys)
% Used with FitBMPSpread
% overlay average pSmad1 vs distance for various times. 

fit_exp = 1;
double_bins = 0;  % limited data BMP500_1flowcell average bins pair wise
keys = strrep(keys, '*', '');
keys = strrep(keys, '_', ' ');
nkeys = length(keys);
nbins = length(stats(1).cnts);
if double_bins
    bins2 = floor(nbins/2);
    for i =1:nkeys
        stats(i).std = (stats(i).std).^2;
        for j = 1:bins2
            stats(i).npts(j) = (stats(i).npts(2*j-1) + stats(i).npts(2*j));
            stats(i).cnts(j) = (stats(i).cnts(2*j-1) + stats(i).cnts(2*j));
            stats(i).std(j) = (stats(i).std(2*j-1) + stats(i).std(2*j))/2;
        end
        stats(i).lims(2:2:end) = [];
        stats(i).npts((bins2+1):end) = [];
        stats(i).cnts((bins2+1):end) = [];
        stats(i).std((bins2+1):end) = [];
        stats(i).std = sqrt(stats(i).std);
    end
    nbins = bins2;
end

avr_smad = zeros(nkeys, nbins);
std_smad = zeros(nkeys, nbins);
for i = 1:nkeys
    avr_smad(i, :) = stats(i).cnts ./ stats(i).npts;
    std_smad(i, :) = stats(i).std ./ stats(i).npts;
end
% Width:  425.0961 microns == 1024 pixels from ImageJ. COnvert bin to mu
% units. 
% nbins = nbins/2;  %% hack plot half the range in mu over x-axis
scl = 425.0961/1024;
bins = scl*floor( (stats(1).lims(1:nbins) + stats(1).lims(2:(nbins+1)))/2 );

% collect data to write to file. Columns: bins, [avr,std,fit] for 3 times
mtx2output = zeros(1+3*nkeys, nbins);  % transpose at end.
ptr=1;
mtx2output(ptr,:) = bins;

figure
symb = {'r', 'g', 'b','k'};
for i = 1:nkeys
    hdl(i) = errorbar(bins, avr_smad(i,1:nbins), std_smad(i,1:nbins), ['*',symb{i}]);
    mtx2output(ptr+1,:) = avr_smad(i,1:nbins);
    mtx2output(ptr+2,:) = std_smad(i,1:nbins);
    hold on
    if fit_exp
        cutoff = 1;  %0.5; %only plot avr_smad > this ratio 
        idx = find(avr_smad(i,1:nbins) < cutoff, 1);
        if isempty(idx)
            idx = nbins;
        end
        if idx > 5
            [poly, smad_fit] = lin_fit(bins(1:idx), avr_smad(i,1:idx));
            plot(bins(1:idx),smad_fit(1:idx), symb{i}); 
            fprintf('data= %d symb= %s, exp decay len,rate mu= %6.4f, %6.4f, intercept= %6.4f\n',...
                i, symb{i}, -1/poly(1), poly);
        else
            smad_fit = zeros(1,nbins);  % dummy value for no fit
        end
        mtx2output(ptr+3,1:idx) = smad_fit(1:idx);
    end
    ptr = ptr + 3;
end
hold off
legend(hdl, keys, 'FontSize',18);
xlabel('Distance to closest secreting cell ({\mu}m)', 'FontSize',18);
ylabel('pSMAD1 intensity', 'FontSize',18);

csvwrite('saved_data.csv', mtx2output');
return

function [poly, fit] = lin_fit(bins, smad)
% fit data over middle of range, return fit to smad from bins(1..)
% only fit smad > cutoff
ifirst = 1;     % ignore the first point
poly = polyfit(bins(ifirst:end), log(smad(ifirst:end)), 1);
fit = exp(poly(1)*bins + poly(2));


% if fit_exp  %% semilog plots only
%     cutoff = 0.5; %only plot avr_smad > this ratio
%     figure,
%     for i = 1:nkeys
%         idx = find(avr_smad(i,:) < cutoff, 1);
%         if isempty(idx)
%             idx = nbins;
%         end
%         semilogy(bins(1:idx), avr_smad(i,1:idx), symb{i});
%         hold on
%     end
%     hold off
%     legend(keys, 'FontSize',18);
%     xlabel('Distance to closest secreting cell ({\mu}m)', 'FontSize',16);
%     
%     ylabel('pSMAD1 intensity', 'FontSize',16);
% end