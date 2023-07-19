function [Image,para] = Reconstruction_IR_SMS(para)

fprintf('\n');
disp('$$$ MRI Reconstruction by Ye $$$');
disp('$$$    phye1988@gmail.com    $$$');
fprintf('\n');



%% load data
load([para.dir.load_kSpace_dir,para.dir.load_kSpace_name]);
kSpace = kSpace*10^8;
sx = size(kSpace,1);
nor = size(kSpace,2);
nc = size(kSpace,5);
keyboard
kSpace = reshape(kSpace,[sx,nor,10,nc]);
coil_drop = [8,9,11,12,18,25];
kSpace(:,:,:,coil_drop) = [];
% order = kSpace_info.InversionTimes_ms(:);
% [~,order] = sort(order);
% order = [6,1,7,2,3,4,5,8,9,10];
order = [1,2,3,4,5];
kSpace = kSpace(:,:,order,:);

[kx,ky] = get_k_coor(sx,kSpace_info.angle_mod(:,:,order),0,round(sx/2+1));

phase(1,:,:,1,2) = exp(1i*kSpace_info.phase_mod(:,:,order)*2*pi/3);
phase(1,:,:,1,3) = conj(phase(:,:,:,:,2));
phase(1,:,:,1,1) = 1;

N = NUFFT.init_new(kx,ky,1.5,[6,6]);

Image = NUFFT.NUFFT_adj_new(kSpace.*conj(phase),N);

Data.kSpace = kSpace;
Data.phase_mod = phase;
Data.N = N;
Data.sens_map = get_sens_map(sum(Image,3),'SMS');
Data.first_est = sum(Image.*conj(Data.sens_map),4);


para.Recon.weight_tTV = max(abs(Data.first_est(:)))*0;
para.Recon.weight_sTV = max(abs(Data.first_est(:)))*0.000;
para.Recon.no_comp = size(Image,4);
para.Recon.noi = 20;
para.setting.ifplot = 1;
para.Recon.type = 'NUFFT';
para.Recon.nSMS = 3;
para.Recon.break = 0;

Image_recon = STCR_conjugate_gradient_T1(Data,para);

keyboard


end