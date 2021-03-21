function check_mask(data, mask)
% overlay boundary of mask on data it was derived from
figure, imshow(data, []);
bndry = bwboundaries(mask);
hold on
for i = 1:length(bndry)
    arc1 = bndry{i};
    plot(arc1(:,2), arc1(:,1), 'r', 'LineWidth',1);
end