function [Image,para] = Reconstruction_old_3D(para)

fprintf('\n');
disp('$$$ 3D MRI Reconstruction by Ye $$$');
disp('$$$      phye1988@gmail.com     $$$');
fprintf('\n');

%%%%% load parameters from para
para = prepare_para(para);

%%%%% load k-space data and PCA on coils
if para.setting.PCA == 0
    load([para.dir.load_kSpace_dir para.dir.load_kSpace_name])
    kSpace = Data.Kspace; clear Data;
    kSpace = kSpace*1e8;
else
    [kSpace,Kx,Ky,para] = PCA_3D_kSpace_yt(para);
end

%%%%% pre-interpolation
if length(size(kSpace)) > 4
    [Data,para] = pre_interp_3D(kSpace,Kx,Ky,para);
end

%%%%% pre-process data 
Data.kSpace(isnan(Data.kSpace)) = 0;
% The following operations are making sure the image has correct phase and
% no need to fftshift at each iteration. The k-space has check board phase
% and is fftshift though.
Data.kSpace = fftshift3(Data.kSpace);
Data.mask = logical(abs(Data.kSpace));

Data.kSpace = ifft3(Data.kSpace);
Data.kSpace = fftshift3(Data.kSpace);
Data.kSpace = fft3(Data.kSpace);
Data.kSpace = Data.kSpace .* Data.mask;

%%%%% sensitivity map and first estimation
%[sens_map,para] = get_sens_map_3D(kSpace,para);
im = ifft3(Data.kSpace);
Data.sens_map = get_sens_map(im,'3D');

para.Recon.nor = squeeze(sum(kSpace(144,:,:,1,1)~=0));
para.Recon.sx = size(kSpace,1);
Data.filter = ramp_filter_for_pre_interp_3D(para);
Data.first_est = sum(bsxfun(@times,ifft3(Data.kSpace.*Data.filter),conj(Data.sens_map)),5);

para.Recon.PD_frames = false(1,size(im,4));
para.Recon.PD_frames(1:10) = true;
para.Recon.bin1 = para.Recon.PD_frames;
para.Recon.bin2 = ~para.Recon.PD_frames;
if sum(para.Recon.bin2) > 90
    para.Recon.bin2(101:end) = false;
end
para.Recon.type = '3D';
para.Recon.no_comp = size(im,5);
sx = size(im,1);
para.Recon.kSpace_size = [sx,sx];
Image = Recon_3D_2_injections(Data,para);
Image = cat(4,Image.pd,Image.inj1);
%% save image and parameter
disp('Saving image into Results...')
Image = gather(Image);
save([para.dir.save_recon_img_mat_dir,para.dir.save_recon_img_name],'Image','para','-v7.3');
disp('Reconstruction done');fprintf('\n')