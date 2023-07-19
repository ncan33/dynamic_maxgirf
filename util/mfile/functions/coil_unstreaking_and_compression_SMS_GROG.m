function Data = coil_unstreaking_and_compression_SMS_GROG(kSpace,kSpace_info)
idx = squeeze(sum(sum(sum(sum(kSpace)))) == 0);
kSpace(:,:,:,:,idx) = [];
kSpace = kSpace * 10^8;

kSpace = permute(kSpace,[1,2,4,3,5]);
kSpace = phase_correction_022718(kSpace);

kSpace = permute(kSpace,[1,2,3,5,4]);
%kSpace_raw = kSpace;
%kSpace = permute(kSpace,[1,2,4,5,3]);
[sx,nor,nof,ns,nc] = size(kSpace);
phase_mod = kSpace_info.phase_mod;
angle_mod = kSpace_info.angle_mod;

phase(1,:,1,1,1,1) = exp(-1i*phase_mod*0);
phase(1,:,1,1,1,2) = exp(-1i*phase_mod*2*pi/3);
phase(1,:,1,1,1,3) = exp(-1i*phase_mod*4*pi/3);
phase = reshape(phase,1,nor,nof,1,1,3);

[kx,ky] = get_k_coor(sx,angle_mod.',0,round((sx+1)/2));
kx = reshape(kx,sx,nor,nof);
ky = reshape(ky,sx,nor,nof);

para.Recon.nSMS = 3;
para.core_size = [1,1];
para.over_sampling = 1;
para.Recon.kSpace_size = [sx,sx];
para.Recon.nor = nor;
para.Recon.sx = sx;

threshold = 25;
keyboard
for i=1:3
    [Data{i}.G,Data{i}.kSpace,~] = GROG.GROG_seperate_SMS_GNUFFT(squeeze(kSpace(:,:,:,i,:)),kx,ky,phase, para);
    s = round(Data{i}.G{1}.core_size/2);
    Data{i}.kSpace = Data{i}.kSpace(s+1:s+Data{i}.G{1}.sx_over, s+1:s+Data{i}.G{1}.sx_over,:,:,:,:,:);
end

for i=1:3
Data{i}.kSpace = fftshift2(Data{i}.kSpace);
Data{i}.mask = logical(abs(Data{i}.kSpace(:,:,:,1,:,:,:)));
Data{i}.filter = ramp_filter_for_pre_interp(para);

Data{i}.kSpace = ifft2(Data{i}.kSpace);
Data{i}.kSpace = fftshift2(Data{i}.kSpace);
Data{i}.kSpace = fft2(Data{i}.kSpace);
Data{i}.kSpace = Data{i}.kSpace .* Data{i}.mask;

kSpace_sms(:,:,:,:,1) = sum(Data{i}.kSpace,7);
kSpace_sms(:,:,:,:,2) = Data{i}.kSpace(:,:,:,:,:,:,1) + exp(-1i*2*pi/3)*Data{i}.kSpace(:,:,:,:,:,:,2) + exp(-1i*4*pi/3)*Data{i}.kSpace(:,:,:,:,:,:,3);
kSpace_sms(:,:,:,:,3) = Data{i}.kSpace(:,:,:,:,:,:,1) + exp(-1i*4*pi/3)*Data{i}.kSpace(:,:,:,:,:,:,2) + exp(-1i*2*pi/3)*Data{i}.kSpace(:,:,:,:,:,:,3);

im = ifft2(kSpace_sms.*Data{i}.filter);
im = fftshift2(im);

im_one = squeeze(mean(im(:,:,11:13,:,:),3));
im_ref = squeeze(mean(im(:,:,11:40,:,:),3));

for j=1:3
    im_one_temp = squeeze(im_one(:,:,:,j));
    im_ref_temp = squeeze(im_ref(:,:,:,j));
    %scale_temp = max(max(abs(im_ref_temp)))./max(max(abs(im_one_temp)));
    %scale_temp = nor_ref/nor_one;
    n_temp = abs(im_one_temp - im_ref_temp).^2;
    d_temp = abs(crop_half_FOV(im_ref_temp)).^2;
    score_temp(j,:) = sum(sum(n_temp))./sum(sum(d_temp));
end

score = squeeze(sum(score_temp,1));
%score = score_temp.';
score = score/min(score(:));
score(score<threshold) = 1;
score(score>threshold) = 1./score(score>threshold)*threshold;
score = permute(score,[3,1,2]);

im_pca = sum(im_ref,4).*score;
im_pca = crop_half_FOV(im_pca);
im_pca = reshape(im_pca,sx*sx/4,nc);


coeff = pca(im_pca);
coeff = permute(coeff(:,1:8),[3,4,5,1,2]);
score = permute(score,[4,1,2,3]);
Data{i}.kSpace = Data{i}.kSpace.*coeff.*score;
Data{i}.kSpace = sum(Data{i}.kSpace,4);
Data{i}.kSpace = permute(Data{i}.kSpace,[1,2,3,5,4,6,7]);
end
