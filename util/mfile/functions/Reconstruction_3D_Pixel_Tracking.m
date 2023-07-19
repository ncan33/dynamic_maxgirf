function [Image,para] = Reconstruction_3D_Pixel_Tracking(para)

fprintf('\n');
disp('$$$ MRI Reconstruction by Ye $$$');
disp('$$$    phye1988@gmail.com    $$$');
fprintf('\n');
%% load parameters from para
para = prepare_para(para);

%% load k-space data
load([para.dir.load_kSpace_dir,para.dir.load_kSpace_name])

%% prepare k-space data
[Data,para] = prepare_Data_3D_new(kSpace,kSpace_info,para);

%% get the STCR recon parameters
scale_image = max(abs(Data.first_est(:)));
para.Recon.weight_tTV = scale_image*para.weight_tTV;
para.Recon.weight_sTV = scale_image*para.weight_sTV;
para.Recon.weight_sliceTV = scale_image*para.weight_sliceTV;
clearvars -except Data para

%% reconstructe image
para.Recon.bin1 = para.Recon.PD_frames;
para.Recon.bin2 = ~para.Recon.PD_frames;
Image = Recon_3D_2_injections_pixel_tracking(Data,para);
Image = cat(4,Image.pd,Image.inj1);

%% save image and parameter
disp('Saving image into Results...')
Image = gather(Image);
save([para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name],'Image','para','-v7.3');
disp('Reconstruction done');fprintf('\n')

end