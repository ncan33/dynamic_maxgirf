function mask_out = refine_mask_edge(mask_in)
[sx,sy] = size(mask_in);
s = regionprops(mask_in,'Centroid');
center = s.Centroid;
boundry = bwboundaries(mask_in);
boundry = boundry{1};

CartCoor = boundry - flip(center);
[theta,r] = cart2pol(CartCoor(:,1),CartCoor(:,2));
r_mean = mean(r);
r(r<r_mean) = r_mean;
%r = r+2;
r = cconv(r,[1/3 1/3 1/3],length(r));
[CartCoor(:,1),CartCoor(:,2)] = pol2cart(theta,r);
CartCoor = round(CartCoor + flip(center));
mask_out = poly2mask(CartCoor(:,2),CartCoor(:,1),sx,sy);
%figure,imagesc(mask_out)
%hold on
%plot(center(1),center(2),'r*')

%r = mean(sqrt(sum((boundry - center).^2,2)));
%[X,Y] = meshgrid(1:sx,1:sy);
%X = X - center(1);
%Y = Y - center(2);
%mask_out = sqrt(X.^2 + Y.^2) <= r;
%mask_out = (mask_out - mask_in) > 0 | mask_in;