function mask_out = fit_mask_to_round(mask_in)
[sx,sy] = size(mask_in);
s = regionprops(mask_in,'Centroid');
center = s.Centroid;
boundry = bwboundaries(mask_in);
boundry = boundry{1};
r = mean(sqrt(sum((boundry - flip(center)).^2,2)));
[X,Y] = meshgrid(1:sx,1:sy);
X = X - center(1);
Y = Y - center(2);
keyboard
mask_out = sqrt(X.^2 + Y.^2) <= r;
mask_out = (mask_out - mask_in) > 0 | mask_in;