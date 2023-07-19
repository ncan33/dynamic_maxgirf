function [sys,dia] = auto_gating_long_axis(Image)
Image = squeeze(abs(crop_half_FOV(Image)));
[sx,sy,nof] = size(Image);

LV = FindLVLongAxis_yt(Image,0);
mask = zeros(sx,sy);
mask(LV(1),LV(2)) = true;
mask = bwdist(mask)<10;
MaxImage = max(Image,[],3);
MaxImage = MaxImage/max(MaxImage(:));
mask = seg_yt(MaxImage,mask,100,1,0,1);
[sys,dia] = get_sys_dia(Image,mask);