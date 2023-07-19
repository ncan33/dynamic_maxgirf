function out = NRMSE(im1,im2)

im1 = abs(im1);
im2 = abs(im2);

d = (im1 - im2).^2;
m1 = mean(im1(:));
out = sqrt(mean(d(:)))/m1;