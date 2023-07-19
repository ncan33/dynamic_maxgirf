function mask_out = refine_mask_edge_new(mask_in)
[sx,sy] = size(mask_in);
s = regionprops(mask_in,'Centroid');
center = s.Centroid;
boundry = bwboundaries(mask_in);
boundry = boundry{1};

CartCoor = boundry - flip(center);

[theta,r] = cart2pol(CartCoor(:,1),CartCoor(:,2));
[pks,locs] = findpeaks(r);
p1 = pks(1);
shift = locs(1);
r = circshift(r,-(shift-1));
r(end+1) = r(1);
[pks,locs] = findpeaks(r);
pks(end+1:end+2) = p1;
locs(end+1) = length(r);
locs(end+1) = 1;
locs = circshift(locs,1);
pks = circshift(pks,1);
np = length(pks);
for i=2:np
    d = (pks(i) - pks(i-1))/(locs(i)-locs(i-1));
    temp = pks(i-1)+d:d:pks(i)-d;
    r(locs(i-1)+1:locs(i)-1) = temp;
end
r(end) = [];
r = circshift(r,shift-1);

[CartCoor(:,1),CartCoor(:,2)] = pol2cart(theta,r);
CartCoor = round(CartCoor + flip(center));
mask_out = poly2mask(CartCoor(:,2),CartCoor(:,1),sx,sy);