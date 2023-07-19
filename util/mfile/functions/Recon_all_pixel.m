function [Image_sys,Image_dia] = Recon_PT_seperate_cardiac_bins(Image_ref,Image_coil,Data,para)


cardiac_signal = self_gating_image_space_4(crop_half_FOV(Image_ref));
respiration_signal = get_resporation_bins(crop_half_FOV(Image_ref));


for i=1:size(Image_ref,4)
    sys = cardiac_signal(:,i) == 1;
    respiration_signal_temp = respiration_signal(sys,i);
    Data_sys.filter = Data{i}.filter;
    Data_sys.kSpace = Data{i}.kSpace(:,:,sys,:,:,:,:);
    Data_sys.mask = Data{i}.mask(:,:,sys,:,:,:,:);
    Data_sys.sens_map = Data{i}.sens_map;
    Data_sys.first_est = Data{i}.first_est(:,:,sys,:,:);
    Data_sys.first_guess = Data{i}.first_guess(:,:,sys,:,:);
    Data_sys.sTV_mask = Data{i}.sTV_mask;
    
    para_sys = para{i};
    para_sys.Recon.bins(1,:) = respiration_signal_temp == 1;
    para_sys.Recon.bins(2,:) = respiration_signal_temp == 2;
    para_sys.Recon.bins(3,:) = respiration_signal_temp == 3;
    para_sys.Recon.bins(4,:) = respiration_signal_temp == 4;
    
    idx = sum(para_sys.Recon.bins,2)<=2;
    para{i}.Recon.bins(idx,:) = [];
        
    para{i}.respiration_signal = respiration_signal;
    para{i}.cardiac_signal = cardiac_signal;
    
    Data_sys.Motion = get_motion_SMS(Image_ref(:,:,sys,i,:)./Image_coil(:,:,:,i,:),1,20);
    Data_sys.Motion_bins = get_motion_SMS_bins(Image_ref(:,:,sys,i,:)./Image_coil(:,:,:,i,:),para_sys,1,20);
    
    para_sys.Recon.weight_tTV = para{i}.Recon.weight_tTV*5;
    Image_sys(:,:,:,i,:) = STCR_conjugate_gradient_pixel_bins(Data_sys,para_sys);
end
