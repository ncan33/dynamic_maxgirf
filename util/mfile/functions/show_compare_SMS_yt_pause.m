function show_compare_SMS_yt_pause(img1,img2)

[sx,sy,nof,nSMS,ns] = size(img1);
img1 = reshape(img1,[sx sy nof nSMS*ns]);
img2 = reshape(img2,[sx sy nof nSMS*ns]);
img1 = permute(img1,[1 2 4 3]);
img2 = permute(img2,[1 2 4 3]);
img1 = reshape(img1,[sx,sy*nSMS*ns,nof]);
img2 = reshape(img2,[sx,sy*nSMS*ns,nof]);
figure
for i=1:size(img1,3)
    im1_temp = img1(:,:,i);
    im2_temp = img2(:,:,i);
    temp_scale = abs(im1_temp(:)).' * pinv(abs(im2_temp(:))).';
    imagesc(abs([im1_temp; temp_scale*im2_temp]))
    colormap gray
    brighten(0.4)
    axis equal
    title(i)

    pause
end
end