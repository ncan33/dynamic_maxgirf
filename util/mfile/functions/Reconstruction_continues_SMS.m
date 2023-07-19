function [Image,para] = Reconstruction_continues_SMS(para)

fprintf('\n');
disp('$$$ MRI Reconstruction by Ye $$$');
disp('$$$    phye1988@gmail.com    $$$');
fprintf('\n');
%% load parameters from para
para = prepare_para(para);

%% load k-space data
load([para.dir.load_kSpace_dir,para.dir.load_kSpace_name])
kSpace = phase_correction_cSMS_041318(kSpace,kSpace_info);
kSpace = yExportKspace(kSpace,kSpace_info.angle_mod,0.5);


%% prepare k-space data

[Data,para] = prepare_Data_cSMS(kSpace,kSpace_info,para);
clearvars -except Data para

%% reconstructe image
[Image_pd,para] = STCR_PD(Data,para);
[Image,para] = STCR_conjugate_gradient(Data,para);

%[Image,para] = STCR_gradient_descent(Data,para);
Image = cat(3,Image_pd,Image);
%Image = Image_pd;
%% save image and parameter
disp('Saving image into Results...')
if para.Recon.crop_half_FOV == 1
    Image = abs(crop_half_FOV(Image));
else
    Image = abs(Image);
end

   
save([para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name],'Image','para','-v7.3');
disp('Reconstruction done');fprintf('\n')

end