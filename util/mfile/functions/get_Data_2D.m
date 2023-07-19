function [Data,para] = get_Data_2D(Data,para)

Data.kSpace = fftshift2(Data.kSpace);
Data.mask = logical(abs(Data.kSpace(:,:,:,1,:,:,:)));
Data.kSpace = ifft2(Data.kSpace);
Data.kSpace = fftshift2(Data.kSpace);
Data.kSpace = fft2(Data.kSpace);
Data.kSpace = Data.kSpace.*Data.mask;
Data.first_est = ifft2(Data.kSpace);
if ~isfield(Data,'sens_map')
    Data.sens_map = get_sens_map(Data.first_est,'2D');
end
para.Recon.kSpace_size = [size(Data.kSpace,1),size(Data.kSpace,2)];
para.Recon.image_size = [size(Data.kSpace,1),size(Data.kSpace,2)];
para.Recon.sx = size(Data.kSpace,1);
Data.filter = ramp_filter_for_pre_interp(para);
Data.first_est = ifft2(Data.kSpace.*Data.filter);
Data.first_est = Data.first_est.*Data.sens_map;
Data.first_est = sum(Data.first_est,4);
