function [sys,dia] = auto_gating_short_axis(Image)
Image = squeeze(abs(crop_half_FOV(Image)));
[sx,sy,nof] = size(Image);

[LV,RV] = FindLVRV_yt(Image,0);
mask = zeros(sx,sy);
mask(RV(1),RV(2)) = true;
mask = bwdist(mask)<10;
MaxImage = max(Image,[],3);
MaxImage = MaxImage/max(MaxImage(:));
mask = seg_yt(MaxImage,mask,100,1,0,1);
[sys,dia] = get_sys_dia(Image,mask);