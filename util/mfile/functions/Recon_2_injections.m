function Image = Recon_2_injections(Data,para)

bin1 = para.Recon.bin1; % pd
bin2 = para.Recon.bin2; % 1 injection
bin3 = para.Recon.bin3; % 2 injection

if isfield(Data,'filter')
    Data_temp.filter = Data.filter;
end
Data_temp.sens_map = Data.sens_map;

Data_temp.kSpace = Data.kSpace(:,:,bin1,:);
Data_temp.first_est = Data.first_est(:,:,bin1);
if isfield(Data,'mask')
    Data_temp.mask = Data.mask(:,:,bin1);
end
scale_image = mean(abs(Data_temp.first_est(:)));
para.Recon.weight_tTV = scale_image*para.weight_tTV;
para.Recon.weight_sTV = scale_image*para.weight_sTV;
Image.pd = STCR_conjugate_gradient(Data_temp,para);
Image.pd = crop_half_FOV(abs(Image.pd));

Data_temp.kSpace = Data.kSpace(:,:,bin2,:);
Data_temp.first_est = Data.first_est(:,:,bin2);
if isfield(Data,'mask')
    Data_temp.mask = Data.mask(:,:,bin2);
end
scale_image = mean(abs(Data_temp.first_est(:)));
para.Recon.weight_tTV = scale_image*para.weight_tTV;
para.Recon.weight_sTV = scale_image*para.weight_sTV;
Image.inj1 = STCR_conjugate_gradient(Data_temp,para);
Image.inj1 = crop_half_FOV(abs(Image.inj1));

Data_temp.kSpace = Data.kSpace(:,:,bin3,:);
Data_temp.first_est = Data.first_est(:,:,bin3);
if isfield(Data,'mask')
    Data_temp.mask = Data.mask(:,:,bin3);
end
scale_image = mean(abs(Data_temp.first_est(:)));
para.Recon.weight_tTV = scale_image*para.weight_tTV;
para.Recon.weight_sTV = scale_image*para.weight_sTV;
Image.inj2 = STCR_conjugate_gradient(Data_temp,para);
Image.inj2 = crop_half_FOV(abs(Image.inj2));
