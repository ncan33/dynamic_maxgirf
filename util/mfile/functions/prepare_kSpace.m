function [kSpace,phase_mod,para] = prepare_kSpace(kSpace_all,theta_all,phase_mod_all,RayPosition,para)
% [kSpace,phase_mod,para] = prepare_kSpace(kSpace_all,theta_all,phase_mod_all,RayPosition,para)
t1 = tic;
[sx,~,no_comp,~,ns] = size(kSpace_all);
nof                 = para.Recon.nof;
nSMS                = para.Recon.nSMS;
selected_rays_start = para.selected_rays_start;
selected_rays_end   = para.selected_rays_end;
ifNUFFT             = para.ifNUFFT;
kCenter             = para.kSpace_center;
interp_method       = para.interp_method;

nor = selected_rays_end - selected_rays_start + 1; para.Recon.nor = nor;
%{
temp = find(RayPosition==1);
RayPosition(temp(end):end) = 0; % delete the last few rays
RayPosition(temp(1):temp(2)-1) = 0; %deleta first few rays

index_1 = find(RayPosition==1);
nor_all = diff(index_1);
nof_low_ray_number = find(nor_all<210);

for i=1:length(nof_low_ray_number)
    RayPosition(index_1(nof_low_ray_number(i)):index_1(nof_low_ray_number(i)+1)-1) = 0;
end
%}
RayIndex = logical((RayPosition > selected_rays_start-1) .* (RayPosition < selected_rays_end+1));
theta = theta_all(RayIndex); nof = length(theta)/nor;

theta = reshape(theta,[1,nor,nof]);
kSpace = kSpace_all(:,RayIndex,:,:,:);
kSpace = reshape(kSpace,[sx,nor,nof,no_comp,1,ns]);
phase_mod = phase_mod_all(1,RayIndex,1,:);
phase_mod = reshape(phase_mod,[1,nor,nof,1,nSMS]);

%%%%% pre-interpolation
disp('Pre-interpolate into Cartesian space...')

matObj = matfile([para.dir.load_kSpace_dir,'/',para.dir.load_kSpace_name]);
varlist = who(matObj,'Kx');
if ~isempty(varlist)
    load([para.dir.load_kSpace_dir,'/',para.dir.load_kSpace_name],'Kx','Ky')
    x_coor = double(squeeze(Kx)); y_coor = double(squeeze(Ky));
else
    [x_coor, y_coor] = get_k_coor(sx,theta,ifNUFFT,kCenter);
end

switch interp_method
    case 'grid3'
        kSpace = pre_interp_radial_to_cart(kSpace,x_coor,y_coor);
        %kSpace = pre_interp_radial_to_cart_oversample(kSpace,x_coor,y_coor,para.over_sampling_factor);
    case'nn'
        kSpace = permute(kSpace,[1 2 4 3 5 6]);
        kSpace_c = zeros([sx sx no_comp nof 1 ns]);
        parfor i=1:ns
            kSpace_c(:,:,:,:,1,i) = interp_nn(kSpace(:,:,:,:,1,i),x_coor, y_coor);
        end
        kSpace = permute(kSpace_c,[1 2 4 3 5 6]); clear kSpace_c
    case 'GROG'
        %kSpace = bsxfun(@times,kSpace,conj(phase_mod));

        kSpace = GROG.GROG(squeeze(kSpace),x_coor,y_coor, 0, 1);
        if nSMS == 1
            kSpace = permute(kSpace,[1 2 3 4 6 5]);
        end
end

%fftshift_mask = single(ones(sx));
%fftshift_mask(1:2:sx,:) = -fftshift_mask(1:2:sx,:);
%fftshift_mask(:,1:2:sx) = -fftshift_mask(:,1:2:sx);

%kSpace = bsxfun(@times,kSpace,fftshift_mask); % fftshift on first 2 dimension can put the result image into right phase. Or the recon-image has to be multiplied by fftshift_mask to be in the right phase. Not doing so makes the image sharper, and the k-space has center in the center.

mask = logical(abs(kSpace));
mask = fftshift(mask,1);
mask = fftshift(mask,2);
kSpace = fftshift(kSpace);
kSpace = ifft2(kSpace);
kSpace = fftshift(kSpace);
kSpace = fft2(kSpace);
kSpace = kSpace.*mask;

phase_mod_all = repmat(phase_mod,[sx,1,1,1,1]); clear phase_mod
phase_mod = single(logical(abs(kSpace(:,:,:,1,1,1))));

if nSMS>1
    phase_mod(:,:,:,:,2:nSMS) = pre_interp_radial_to_cart(phase_mod_all(:,:,:,:,2:nSMS),x_coor,y_coor);
    phase_mod(:,:,:,:,2:nSMS) = fftshift(phase_mod(:,:,:,:,2:nSMS),1);
    phase_mod(:,:,:,:,2:nSMS) = fftshift(phase_mod(:,:,:,:,2:nSMS),2);
end

if para.image_orintation == 0
    orin = orintation_detection(abs(fftshift(ifft2(sum(sum(kSpace,3),4)))));
    kSpace = orintate_image(kSpace,orin);
else
    kSpace = orintate_image(kSpace,para.image_orintation);
    phase_mod = orintate_image(phase_mod,para.image_orintation);
end
    
%kSpace = fftshift(fftshift(kSpace,1),2);
%phase_mod = fftshift(fftshift(phase_mod,1),2);
para.CPUtime.prepare_kSpace = toc(t1);toc(t1);fprintf('\n');
  