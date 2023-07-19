function [Image,para] = Reconstruction_ktGROG(para)

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
keyboard
%% reconstructe image
tic
k = ktGROG.Fill_kSpace_KWIC(Data.kSpace);
Image = ifft2(k).*conj(Data.sens_map);
Image = sum(Image,4);
Image = crop_half_FOV(Image);

%Image = tv_denoise(Image,mean(abs(Image(:)))*0.01,50);
for i=1:size(Image,3)
    im_r(:,:,i) = medfilt2(real(Image(:,:,i)),[2,2]);
    im_i(:,:,i) = medfilt2(imag(Image(:,:,i)),[2,2]);
end
Image = sqrt(im_r.^2+im_i.^2);
Image = gather(Image);
toc
%Image = Image_pd;
%% save image and parameter
disp('Saving image into Results...')
%if para.Recon.crop_half_FOV == 1
%    Image = abs(crop_half_FOV(Image));
%else
    Image = abs(Image);
%end

save([para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name],'Image','para','-v7.3');
disp('Reconstruction done');fprintf('\n')

end