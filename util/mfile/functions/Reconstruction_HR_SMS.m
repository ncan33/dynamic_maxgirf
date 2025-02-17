function [Image,para] = Reconstruction_HR_SMS(para)

fprintf('\n');
disp('$$$ MRI Reconstruction by Ye $$$');
disp('$$$    phye1988@gmail.com    $$$');
fprintf('\n');

%% load parameters from para
para_temp = prepare_para(para);clear para

%% load k-space data
[kSpace_all,RayPosition,para_temp] = PCA_kSpace_yt(para_temp);

for i=1:size(kSpace_all,5)
    para{i} = para_temp;
    [Data{i},para{i}] = prepare_Data(kSpace_all(:,:,:,:,i),RayPosition,para{i});
    [Data{i},para{i}] = remove_PD(Data{i},para{i});
    scale_image = mean(abs(Data{i}.first_est(:)));
    para{i}.Recon.weight_tTV = scale_image*para{i}.weight_tTV;
    para{i}.Recon.weight_sTV = scale_image*para{i}.weight_sTV;
end

%% reconstructe PD
for i=1:length(Data)

    noi_temp = para{i}.Recon.noi;

    para{i}.Recon.noi = 50;
    %para{i}.Recon.weight_tTV = para{i}.Recon.weight_tTV/3;
    Image_ref(:,:,:,i,:) = STCR_conjugate_gradient(Data{i},para{i});
    Data{i}.Motion = get_motion_SMS(Image_ref(:,:,:,i,:),para{i}.Recon.motion_noi);
    %para{i}.Recon.weight_tTV = para{i}.Recon.weight_tTV*3;
    para{i}.Recon.noi = noi_temp;
    Data{i}.first_guess = Image_ref(:,:,:,i,:);
    [Image(:,:,:,i,:),para{i}] = STCR_conjugate_gradient_pixel(Data{i},para{i});
    Data{i} = rmfield(Data{i},'first_guess');

end
keyboard
%dia = cardiac_signal<mean(cardiac_signal);
%rep = respiration_signal>2500;

para{1}.Recon.bin1 = dia&rep;
para{1}.Recon.bin2 = dia&~rep;
para{1}.Recon.bin3 = ~dia&rep;
para{1}.Recon.bin4 = ~dia&~rep;

Image_4_bins = Recon_4_bins_SMS_pixel(Data{i},para{i});


%% save image and parameter

if para{1}.Recon.crop_half_FOV == 1
    Image = abs(crop_half_FOV(Image));
    Image_ref = abs(crop_half_FOV(Image_ref));
else
    Image = abs(Image);
    Image_ref = abs(Image_ref);
end

disp('Saving image into Results...')

Image = gather(Image);
Image_ref = gather(Image_ref);
save([para{1}.dir.save_recon_img_mat_dir,para{1}.dir.save_recon_img_name],'Image','Image_ref','para','-v7.3');

disp('Reconstruction done');fprintf('\n')

end