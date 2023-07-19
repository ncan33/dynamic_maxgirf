function [Image,para] = Reconstruction_3D_ungated_102218_phantom(para)

fprintf('\n');
disp('$$$ MRI Reconstruction by Ye $$$');
disp('$$$    phye1988@gmail.com    $$$');
fprintf('\n');
%% load parameters from para
para = prepare_para(para);

%% load k-space data
load([para.dir.load_kSpace_dir,para.dir.load_kSpace_name])
%% prepare k-space data
kSpace = yExportKspace(kSpace,0,0.5);
keyboard
[Data,para] = prepare_Data_3D_low_res_ungated_phantom(kSpace,kSpace_info,para);
%[Data,para] = prepare_Data_3D_ungated_no_self_gating(kSpace,kSpace_info,para);
Image_PD = get_PD_3D_continues_multiple_frame(kSpace,kSpace_info,para);
%[Data,para] = prepare_Data_3D_ungated(Data,kSpace,kSpace_info,para);
% keyboard
%% get the STCR recon parameters
% scale_image = max(abs(Data.first_est(:)));
% para.Recon.weight_tTV = scale_image*para.weight_tTV;
% para.Recon.weight_sTV = scale_image*para.weight_sTV;
% para.Recon.weight_sliceTV = scale_image*para.weight_sliceTV;
% clearvars -except Data para Image_PD

%% reconstructe image
Image = STCR_conjugate_gradient_3D(Data,para);
Image = crop_half_FOV(abs(Image));
% nof = 100;
% nsection = ceil(size(Data.kSpace,4)/nof);
% for i=1:nsection
%     if i~=nsection
%         frames = (1:nof) + (i-1)*nof;
%     else
%         frames = (i-1)*nof:size(kSpace,4);
%     end
%     Data_temp = Data;
%     Data_temp.kSpace = Data_temp.kSpace(:,:,:,frames,:);
%     Data_temp.first_est = Data_temp.first_est(:,:,:,frames);
%     Data_temp.mask = Data_temp.mask(:,:,:,frames);
%     scale_image = mean(abs(Data_temp.first_est(:)));
%     para.Recon.weight_tTV = scale_image*para.weight_tTV;
%     para.Recon.weight_sTV = scale_image*para.weight_sTV;
%     para.Recon.weight_sliceTV = scale_image*para.weight_sliceTV;
%     Image_temp = STCR_conjugate_gradient_3D(Data_temp,para);
%     Image_temp = gather(crop_half_FOV(abs(Image_temp)));
%     if i==1
%         Image = Image_temp;
%     else
%         Image = cat(4,Image,Image_temp);
%     end
% end

%% save image and parameter
disp('Saving image into Results...')
save([para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name],'Image','Image_PD','para','-v7.3');
disp('Reconstruction done');fprintf('\n')

end