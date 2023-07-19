function [im_reg,shifts] = rigid_reg(Image,crop_p,ref_frame)

im_small = Image(crop_p(1):crop_p(2),crop_p(3):crop_p(4),:);
[sx,sy,nof] = size(im_small);
im_reg_small = zeros(size(im_small),'like',im_small);
im_reg = zeros(size(Image),'like',im_small);
shifts = zeros(2,nof);

im_reg_small(:,:,ref_frame) = im_small(:,:,ref_frame);
im_reg(:,:,ref_frame) = Image(:,:,ref_frame);
if ref_frame ~= 1
    for i = ref_frame-1 : -1 : 1
        im_fixed = im_reg_small(:,:,i+1);
        im_moving = im_small(:,:,i);
        xcorr_temp = xcorr2(im_fixed,im_moving);
        [x,y] = find(xcorr_temp == max(max(xcorr_temp)));
        shifts(:,i) = [x-sx,y-sy];
        im_reg(:,:,i) = circshift(Image(:,:,i),shifts(:,i));
        im_reg_small(:,:,i) = im_reg(crop_p(1):crop_p(2),crop_p(3):crop_p(4),i);
    end
end

if ref_frame ~= nof
    for i = ref_frame+1 : nof
        im_fixed = im_reg_small(:,:,i-1);
        im_moving = im_small(:,:,i);
        xcorr_temp = xcorr2(im_fixed,im_moving);
        [x,y] = find(xcorr_temp == max(max(xcorr_temp)));
        shifts(:,i) = [x-sx,y-sy];
        im_reg(:,:,i) = circshift(Image(:,:,i),shifts(:,i));
        im_reg_small(:,:,i) = im_reg(crop_p(1):crop_p(2),crop_p(3):crop_p(4),i);
    end
end