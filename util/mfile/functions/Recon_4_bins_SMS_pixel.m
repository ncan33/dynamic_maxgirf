function Image = Recon_4_bins_SMS_pixel(Data,para)

bin1 = para.Recon.bin1;
bin2 = para.Recon.bin2;
bin3 = para.Recon.bin3;
bin4 = para.Recon.bin4;

Data_temp.filter = Data.filter;
Data_temp.sens_map = Data.sens_map;

Data_temp.kSpace = Data.kSpace(:,:,bin1,:,:,:,:);
Data_temp.first_est = Data.first_est(:,:,bin1,:,:);
Data_temp.mask = Data.mask(:,:,bin1,:,:,:,:);
scale_image = mean(abs(Data_temp.first_est(:)));
para.Recon.weight_tTV = scale_image*para.weight_tTV;
para.Recon.weight_sTV = scale_image*para.weight_sTV;
Image_ref = STCR_conjugate_gradient(Data_temp,para);
Data_temp.Motion = get_motion_SMS(Image_ref,para.Recon.motion_noi);
Image.bin1 = STCR_conjugate_gradient_pixel(Data_temp,para);


Data_temp.kSpace = Data.kSpace(:,:,bin2,:,:,:,:);
Data_temp.first_est = Data.first_est(:,:,bin2,:,:);
Data_temp.mask = Data.mask(:,:,bin2,:,:,:,:);
scale_image = mean(abs(Data_temp.first_est(:)));
para.Recon.weight_tTV = scale_image*para.weight_tTV;
para.Recon.weight_sTV = scale_image*para.weight_sTV;
Image_ref = STCR_conjugate_gradient(Data_temp,para);
Data_temp.Motion = get_motion_SMS(Image_ref,para.Recon.motion_noi);
Image.bin2 = STCR_conjugate_gradient_pixel(Data_temp,para);


Data_temp.kSpace = Data.kSpace(:,:,bin3,:,:,:,:);
Data_temp.first_est = Data.first_est(:,:,bin3,:,:);
Data_temp.mask = Data.mask(:,:,bin3,:,:,:,:);
scale_image = mean(abs(Data_temp.first_est(:)));
para.Recon.weight_tTV = scale_image*para.weight_tTV;
para.Recon.weight_sTV = scale_image*para.weight_sTV;
Image_ref = STCR_conjugate_gradient(Data_temp,para);
Data_temp.Motion = get_motion_SMS(Image_ref,para.Recon.motion_noi);
Image.bin3 = STCR_conjugate_gradient_pixel(Data_temp,para);


Data_temp.kSpace = Data.kSpace(:,:,bin4,:,:,:,:);
Data_temp.first_est = Data.first_est(:,:,bin4,:,:);
Data_temp.mask = Data.mask(:,:,bin4,:,:,:,:);
scale_image = mean(abs(Data_temp.first_est(:)));
para.Recon.weight_tTV = scale_image*para.weight_tTV;
para.Recon.weight_sTV = scale_image*para.weight_sTV;
Image_ref = STCR_conjugate_gradient(Data_temp,para);
Data_temp.Motion = get_motion_SMS(Image_ref,para.Recon.motion_noi);
Image.bin4 = STCR_conjugate_gradient_pixel(Data_temp,para);