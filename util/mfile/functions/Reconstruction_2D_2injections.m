function [Image,para] = Reconstruction_2D_2injections(para)

fprintf('\n');
disp('$$$ MRI Reconstruction by Ye $$$');
disp('$$$    phye1988@gmail.com    $$$');
fprintf('\n');
%% load parameters from para
para = prepare_para(para);

%% load k-space data
[kSpace_all,RayPosition,para] = PCA_kSpace_yt(para);

%% prepare k-space data
[Data,para] = prepare_Data(kSpace_all,RayPosition,para);
clearvars -except Data para

%% reconstructe image
Image = Recon_2_injections(Data,para);

%% save image and parameter
disp('Saving image into Results...')

init = abs(crop_half_FOV(Data.first_est));

save([para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name],'Image','para','init','-v7.3');
disp('Reconstruction done');fprintf('\n')

end