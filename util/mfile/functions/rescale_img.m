function im2 = rescale_img(im1,im2)

[sx,sy,nof,nSMS,ns] = size(im1);
im1 = reshape(im1,[sx sy nof*nSMS*ns]);
im2 = reshape(im2,[sx sy nof*nSMS*ns]);
im2_out = zeros([sx sy nof*nSMS*ns]);
if isa(im1,'gpuArray')
    im2_out = gpuArray(im2_out);
end
for i=1:nof*nSMS*ns
    im1_temp = im1(:,:,i);
    im2_temp = im2(:,:,i);
    temp_scale = abs(im1_temp(:)).' * pinv(abs(im2_temp(:))).';
    im2_out(:,:,i) = im2_temp*temp_scale;
end

im2 = reshape(im2_out,[sx,sy,nof,nSMS,ns]);