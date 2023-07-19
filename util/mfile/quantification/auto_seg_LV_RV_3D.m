function [Image_reg,Image_PD,mask_LV_all]=auto_seg_LV_RV_3D(Image,Image_PD)
keyboard
Image = squeeze(abs(Image));
[sx,sy,nslice,nof] = size(Image);

mid_slice = round(nslice/2);

[LV_mid,RV_mid] = find_LV(Image(:,:,mid_slice,:));

%[mask_LV] = seg_LV_RV_Model_Based(Image(:,:,mid_slice,:),LV_mid);
mask_LV =cRV(Image(:,:,mid_slice,:),LV_mid);
%mask_LV_all = seg_Ventricle_3D_model_based(Image,mid_slice,mask_LV);
mask_LV_all = seg_Ventricle_3D(Image,mid_slice,mask_LV);

% for i=1:5
%     Image_reg = reg_Ventricle_3D_model_based(Image,mask_LV_all);
%     mask_LV_all = seg_Ventricle_3D_MB_step_2(Image_reg,mask_LV_all);
% end

for i=1:5
    Image_reg = reg_Ventricle_3D(Image,mask_LV_all);
    mask_LV_all = seg_Ventricle_3D_step_2(Image_reg,mask_LV_all);
end

mask_larger = bwdist(mask_LV_all)<20;
if exist('Image_PD')
    Image_PD = rigid_reg_PD_3D_yt(Image_PD,Image,mask_larger);
    Image_PD(Image_PD<1) = 1;
else
    Image_PD = [];
end
return

for i=1:5
    mask_Myo_all = seg_Myo_3D_MB(Image_reg,mask_LV_all,0.05);
    Image_reg = reg_Myo_3D_MB(Image,mask_LV_all,mask_Myo_all);
end

% Image_PD = rigid_reg_PD_3D_yt(Image_PD,Image,mask_Myo_all|mask_LV_all);
% Image_PD(Image_PD<1) = 1;
% return


flip_angle = 15;
TR = 6;
PD_flip_angle = 2;
Number_of_PD_rays = 1800;

SI_curves = Image_reg./Image_PD;
SI_curves = reshape(SI_curves,sx*sy*nslice,nof);

GD_curves = Quant.SI2Gd_SPGR(SI_curves',flip_angle,TR,PD_flip_angle,Number_of_PD_rays);
GD_curves = reshape(GD_curves',[sx,sy,nslice,nof]);
GD_curves = permute(GD_curves,[4,1,2,3]);


TimeStamp = para.kSpace_info.TimeStamp(para.Recon.bins(1,:));

GD_curves_second = interp1(TimeStamp,GD_curves,ceil(min(TimeStamp)):floor(max(TimeStamp)));
GD_curves_second = permute(GD_curves_second,[2,3,4,1]);


QuantMap = Blind_Quant(GD_curves_second,mask_LV_all,mask_Myo_all);