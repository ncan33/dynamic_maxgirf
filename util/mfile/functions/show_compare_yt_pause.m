function show_compare_yt_pause(img1,img2,slice)
if nargin<3
    slice = 1;
end
img1 = squeeze(img1);
img2 = squeeze(img2);

[sx,sy,nof,nSMS,ns] = size(img1);
img1 = reshape(img1,[sx sy nof nSMS*ns]);
img2 = reshape(img2,[sx sy nof nSMS*ns]);
figure
for i=1:size(img1,3)
    im1_temp = img1(:,:,i,slice);
    im2_temp = img2(:,:,i,slice);
    temp_scale = abs(im1_temp(:)).' * pinv(abs(im2_temp(:))).';
    imagesc(abs([im1_temp temp_scale*im2_temp]))
    colormap gray
    brighten(0.4)
    axis equal
    title(i)

    pause
end
end