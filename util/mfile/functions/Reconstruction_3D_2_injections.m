function [Image,para] = Reconstruction_3D_2_injections(para)

fprintf('\n');
disp('$$$ MRI Reconstruction by Ye $$$');
disp('$$$    phye1988@gmail.com    $$$');
fprintf('\n');
%% load parameters from para
para = prepare_para(para);

%% load k-space data
[kSpace,kx,ky,para] = PCA_3D_kSpace_yt(para);

%% prepare k-space data
[Data,para] = prepare_Data_3D(kSpace,kx,ky,para);
clearvars -except Data para

%% reconstructe image
Image = Recon_3D_2_injections(Data,para);

%% save image and parameter
disp('Saving image into Results...')

save([para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name],'Image','para','-v7.3');
disp('Reconstruction done');fprintf('\n')

end