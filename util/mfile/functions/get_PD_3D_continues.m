function Image_PD = get_PD_3D_continues(kSpace_radial,kSpace_info,para)

kSpace_radial = kSpace_radial(:,1:kSpace_info.NumberOfPDlines,:);
theta_all = kSpace_info.angle_mod(1:kSpace_info.NumberOfPDlines);
slice_idx = kSpace_info.slice_idx(1:kSpace_info.NumberOfPDlines);
if min(slice_idx(:)) <= 0
    slice_idx = slice_idx - min(slice_idx) + 1;
end
sx = size(kSpace_radial,1);
nos = max(slice_idx(:)) + 2;
nc = size(kSpace_radial,3);

kSpace = zeros(sx,1,nos,1,nc);
theta = zeros(1,1,nos);

for islice=1:nos-2
    slice_idx_temp = slice_idx == islice;
    kSpace(:,1:sum(slice_idx_temp),islice,1,:) = kSpace_radial(:,slice_idx_temp,:);
    theta(:,1:sum(slice_idx_temp),islice) = theta_all(:,slice_idx_temp,:);
end

[kx,ky] = get_k_coor(sx,theta,0,round(sx/2)+1);
Data.kSpace = GROG.GROG_3D(kSpace,kx,ky,0,1);
Data.kSpace = fftshift3(Data.kSpace);
Data.mask = logical(abs(Data.kSpace));
Data.first_est = ifft3(Data.kSpace);
Data.first_est = fftshift3(Data.first_est);
Data.kSpace = fft3(Data.first_est);
Data.kSpace = Data.kSpace.*Data.mask;
Data.sens_map = get_sens_map(Data.first_est,'3D');
Data.first_est = Data.first_est.*conj(Data.sens_map);
Data.first_est = sum(Data.first_est,5);
Data.filter = 1;

scale = max(abs(Data.first_est(:)));
para.Recon.weight_tTV = 0;
para.Recon.weight_sTV = scale*para.weight_sTV;
para.Recon.weight_sliceTV = scale*para.weight_sliceTV;
para.Recon.epsilon = eps('single');
para.Recon.step_size = 2;
para.Recon.ifContinue = 0;
para.Recon.no_comp = 8;
para.Recon.noi = 100;
para.Recon.type = '3D';
para.Recon.break = 0;

Image_PD = STCR_conjugate_gradient_3D(Data,para);
Image_PD = abs(crop_half_FOV(Image_PD));