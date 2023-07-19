function [Image_sys,Image_dia,para] = Reconstruction_SMS_multi_slice_ungated_PT_seperate_cardiac_phase(para)

fprintf('\n');
disp('$$$ MRI Reconstruction by Ye $$$');
disp('$$$    phye1988@gmail.com    $$$');
fprintf('\n');

%% load parameters from para
para_temp = prepare_para(para);clear para

%% load k-space data
if para_temp.setting.unstreaking
    [kSpace_all,RayPosition,para_temp] = PCA_kSpace_unstreaking_yt(para_temp);
else
    [kSpace_all,RayPosition,para_temp] = PCA_kSpace_yt(para_temp);
end

for i=1:size(kSpace_all,5)
    para{i} = para_temp;
    if size(para_temp.phase_mod,2) > 1
        para{i}.phase_mod = para{i}.phase_mod(:,i);
        para{i}.angle_mod = para{i}.angle_mod(:,i);
    end
    if length(para_temp.image_orintation)>1
        para{i}.image_orintation = para_temp.image_orintation(i);
    end
    [Data{i},para{i}] = prepare_Data(kSpace_all(:,:,:,:,i),RayPosition,para{i});
end

%% reconstructe PD
for i=1:length(Data)
    [Image_PD_ref(:,:,:,i,:),para{i}] = STCR_PD(Data{i},para{i});
    para{i}.Recon.noi = 40;
    Data{i}.Motion = get_motion_SMS(Image_PD_ref(:,:,:,i,:),1,20);
    Data{i}.first_guess = Image_PD_ref(:,:,:,i,:);
    [Image_PD(:,:,:,i,:),para{i}] = STCR_PD_pixel(Data{i},para{i});
    [Data{i},para{i}] = remove_PD(Data{i},para{i});
end

[Image_coil,sTV_mask] = get_coil_shading_map(max(abs(Image_PD),[],3));

%% reconstruct ref
for i=1:length(Data)
    para{i}.Recon.noi = para_temp.Recon.noi;
    para{i}.step_size = para_temp.step_size;
    para{i}.Recon.weight_tTV = para{i}.Recon.weight_tTV/5;
    [Image_ref(:,:,:,i,:),para{i}] = STCR_conjugate_gradient(Data{i},para{i});
    Data{i}.first_guess = Image_ref(:,:,:,i,:);
    para{i}.Recon.noi = 40;
    if para_temp.setting.adaptive_sTV
        Data{i}.sTV_mask = sTV_mask(:,:,1,i,:);
    else
        Data{i}.sTV_mask = 1;
    end
end

[Image_sys,Image_dia,para] = Recon_PT_seperate_cardiac_bins(Image_ref,Image_coil,Data,para);

%% save image and parameter
%{
if isempty(Image_PD)
elseif isempty(Image)
    Image = Image_PD;
else
    Image = cat(3,Image_PD,Image);
    Image_ref = cat(3,Image_PD_ref,Image_ref);
end
%}
if para{1}.Recon.crop_half_FOV == 1
    Image_coil = crop_half_FOV(Image_coil);
    Image_sys = abs(crop_half_FOV(Image_sys));
    Image_dia = abs(crop_half_FOV(Image_dia));
    Image_ref = abs(crop_half_FOV(Image_ref));
else
    Image_sys = abs(Image_sys);
    Image_dia = abs(Image_dia);
    Image_ref = abs(Image_ref);
end

disp('Saving image into Results...')

Image_sys = gather(Image_sys);
Image_dia = gather(Image_dia);
save([para{1}.dir.save_recon_img_mat_dir,para{1}.dir.save_recon_img_name],'Image_sys','Image_dia','Image_PD','Image_ref','Image_coil','para','-v7.3');

disp('Reconstruction done');fprintf('\n')

end