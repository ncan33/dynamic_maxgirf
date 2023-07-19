function [Image_sys,Image_dia,para] = Reconstruction_SMS_multi_slice_ungated_seperate_bins(para)

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
    if length(para_temp.image_orintation) > 1
        para{i}.image_orintation = para_temp.image_orintation(i);
    end
    [Data{i},para{i}] = prepare_Data(kSpace_all(:,:,:,:,i),RayPosition,para{i});
end

%% reconstructe PD and Ref
%parpool(3)
for i=1:length(Data)
    para{i}.Recon.noi = 100;
    Image_PD(:,:,:,i,:) = STCR_PD(Data{i},para{i});  
    [Data{i},para{i}] = remove_PD(Data{i},para{i});
    para{i}.Recon.noi = 30;
    para{i}.Recon.weight_tTV = para{i}.Recon.weight_tTV/5;
    Image_ref(:,:,:,i,:) = STCR_conjugate_gradient(Data{i},para{i});
    para{i}.Recon.noi = 100;
    para{i}.Recon.weight_tTV = para{i}.Recon.weight_tTV/5;
end

%% bin data
para = auto_gating_cardiac_SMS_yt(Image_ref,para);

%% seperate bins STCR
%Image = zeros(size(Image_ref),'like',Image_ref);
nof = size(Image_ref,3);
for i=1:length(Data)
    nof = min([nof;sum(para{i}.Recon.bins(1:2,:),2)]);
end

for i=1:length(Data)

    Data_temp.filter = Data{i}.filter;
    Data_temp.sens_map = Data{i}.sens_map;
    
    % sys
    bin_temp = find(para{i}.Recon.bins(1,:));
    bin_temp = bin_temp(1:nof);
    Data_temp.kSpace = Data{i}.kSpace(:,:,bin_temp,:,:,:,:);
    Data_temp.first_est = Data{i}.first_est(:,:,bin_temp,:,:);
    Data_temp.mask = Data{i}.mask(:,:,bin_temp,:,:,:,:,:);
    Image_sys(:,:,:,i,:) = STCR_conjugate_gradient(Data_temp,para{i});
    
    % dia
    bin_temp = find(para{i}.Recon.bins(2,:));
    bin_temp = bin_temp(1:nof);
    Data_temp.kSpace = Data{i}.kSpace(:,:,bin_temp,:,:,:,:);
    Data_temp.first_est = Data{i}.first_est(:,:,bin_temp,:,:);
    Data_temp.mask = Data{i}.mask(:,:,bin_temp,:,:,:,:,:);
    Image_dia(:,:,:,i,:) = STCR_conjugate_gradient(Data_temp,para{i});

%     for j=1:length(para{i}.Recon.bins(:,1))
%         bin_temp = para{i}.Recon.bins(j,:);
%         Data_temp.kSpace = Data{i}.kSpace(:,:,bin_temp,:,:,:,:);
%         Data_temp.first_est = Data{i}.first_est(:,:,bin_temp,:,:);
%         Data_temp.mask = Data{i}.mask(:,:,bin_temp,:,:,:,:,:);
%         Image(:,:,bin_temp,i,:) = STCR_conjugate_gradient(Data_temp,para{i});
%     end

end

% if isempty(Image_PD)
% elseif isempty(Image)
%     Image = Image_PD;
% else
%     Image = cat(3,Image_PD,Image);
%     %Image_ref = cat(3,Image_PD_ref,Image_ref);
% end

%% save image and parameter

if para{1}.Recon.crop_half_FOV == 1
    Image_PD = abs(crop_half_FOV(Image_PD));
    Image_sys = abs(crop_half_FOV(Image_sys));
    Image_dia = abs(crop_half_FOV(Image_dia));
    %Image_ref = abs(crop_half_FOV(Image_ref));
else
    Image_PD = abs(Image_PD);
    Image_sys = abs(Image_sys);
    Image_dia = abs(Image_dia);
    %Image_ref = abs(Image_ref);
end

disp('Saving image into Results...')
% 
% Image = gather(Image);
save([para{1}.dir.save_recon_img_mat_dir,para{1}.dir.save_recon_img_name],'Image_PD','Image_sys','Image_dia','para','-v7.3');

disp('Reconstruction done');fprintf('\n')

end