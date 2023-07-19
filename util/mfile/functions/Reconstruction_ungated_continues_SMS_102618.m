function [Image_sys,para] = Reconstruction_ungated_continues_SMS_102618(para)

fprintf('\n');
disp('$$$ MRI Reconstruction by Ye $$$');
disp('$$$    phye1988@gmail.com    $$$');
fprintf('\n');
%% load parameters from para
para = prepare_para(para);

%% load k-space data
load([para.dir.load_kSpace_dir,para.dir.load_kSpace_name])
%kSpace = phase_correction_cSMS_041318(kSpace,kSpace_info);

%% prepare k-space data

%kSpace = yExportKspace(kSpace,1,0.5);
%test_pseudo_SMS_grappa(kSpace,kSpace_info,para);
%recon_2D_ungated_SPGR_multi_set_SMS_low_rank(kSpace,kSpace_info,para);
recon_2D_ungated_SPGR_multi_set_SMS_low_rank_car_phase_resolved(kSpace,kSpace_info,para);
%[Image_sys,Image_dia,Image_cine,Image_ref,para] = prepare_Data_2D_ungated_slice_group_4_MB_4(kSpace,kSpace_info,para);
% clearvars -except Data para

%% reconstructe image
% [Image_pd,para] = STCR_PD(Data,para);
% [Image,para] = STCR_conjugate_gradient(Data,para);

%[Image,para] = STCR_gradient_descent(Data,para);
% Image = cat(3,Image_pd,Image);
%Image = Image_pd;
%% save image and parameter
% disp('Saving image into Results...')
% if para.Recon.crop_half_FOV == 1
%     Image_sys = abs(crop_half_FOV(Image_sys));
%     Image_dia = abs(crop_half_FOV(Image_dia));
%     Image_cine = abs(crop_half_FOV(Image_cine));
%     Image_ref = abs(crop_half_FOV(Image_ref));
% else
%     Image_sys = abs(Image_sys);
%     Image_dia = abs(Image_dia);
%     Image_cine = abs(Image_cine);
%     Image_ref = abs(Image_ref);
% end
% 
%    
% save([para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name],'Image_sys','Image_dia','Image_cine','Image_ref','para','-v7.3');
% disp('Reconstruction done');fprintf('\n')

end