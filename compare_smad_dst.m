function compare_smad_dst()
%
% run in editor mode, to edit files in FitFlowCellData, read in two data sets
% and then compare with call to KS2

scaled_smad = FitFlowCellData(0);
dist1 = dist2compare(scaled_smad);
figure, histogram(dist1, 15);

scaled_smad = FitFlowCellData(0);
dist2 = dist2compare(scaled_smad);
figure, histogram(dist2, 15);

[~, prob, kstest] = kstest2(dist1, dist2);
fprintf('kstest2:  prob= %6.2e, kstest= %6.2e\n', prob, kstest);
return

%%%%%%
function dist = dist2compare(scaled_smad)
% from stats_smad_dst, take scaled_smad(nbins, #noggin cells) and nomalize radial
% distributions relative to inner most ring

mx = max(scaled_smad);
dist = mx - scaled_smad(1,:);