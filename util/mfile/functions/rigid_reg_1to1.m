function im_reg = rigid_reg_1to1(im_moving,im_fixed)

[sx,sy] = size(im_fixed);
xcorr_temp = xcorr2(im_moving,im_fixed);
[x,y] = find(xcorr_temp == max(max(xcorr_temp)));
shifts = [x-sx,y-sy];
im_reg = circshift(im_moving,shifts);