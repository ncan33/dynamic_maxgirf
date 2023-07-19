function [kSpace,para] = prepare_kSpace_sSMS(kSpace_all,theta_all,phase_mod_all,RayPosition,para)
% [kSpace,phase_mod,para] = prepare_kSpace(kSpace_all,theta_all,phase_mod_all,RayPosition,para)
t1 = tic;
[sx,~,no_comp,~,ns] = size(kSpace_all);
%nof                 = para.Recon.nof;
nSMS                = para.Recon.nSMS;
selected_rays_start = para.selected_rays_start;
selected_rays_end   = para.selected_rays_end;
ifNUFFT             = para.ifNUFFT;
kCenter             = para.kSpace_center;
interp_method       = para.interp_method;

nor = selected_rays_end - selected_rays_start + 1; para.Recon.nor = nor;
RayIndex = logical((RayPosition > selected_rays_start-1) .* (RayPosition < selected_rays_end+1));
theta = theta_all(RayIndex); nof = length(theta)/nor;
theta = reshape(theta,[1,nor,nof]);
kSpace = kSpace_all(:,RayIndex,:,:,:);
kSpace = reshape(kSpace,[sx,nor,nof,no_comp,1,ns]);
phase_mod = phase_mod_all(1,RayIndex,1,:);
phase_mod = reshape(phase_mod,[1,nor,nof,1,nSMS]);

%%%%% pre-interpolation
disp('Pre-interpolate into Cartesian space...')

matObj = matfile(para.Recon.PCA_dir);
varlist = who(matObj,'Kx');
if para.asymmetry
    sx = para.AsymmetryRayEnd;
end
if ~isempty(varlist)
    load(para.Recon.PCA_dir,'Kx','Ky')
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

%        [x_all, y_all] = get_k_coor(sx,theta_all,ifNUFFT,kCenter);
%        x_all = reshape(x_all,[288 340 70]); y_all = reshape(y_all,[288 340 70]);
%        kSpace = GROG.GROG_seperate_SMS_train_core_using_all_data(reshape(kSpace_all,[288 340 70 8]),squeeze(kSpace),x_coor,y_coor,x_all,y_all, 0, 1);
        kSpace = GROG.GROG_seperate_SMS(squeeze(kSpace),x_coor,y_coor,phase_mod(:,:,:,:,2), 0, 1, para);
end

fftshift_mask = single(ones(sx));
fftshift_mask(1:2:sx,:) = -fftshift_mask(1:2:sx,:);
fftshift_mask(:,1:2:sx) = -fftshift_mask(:,1:2:sx);

kSpace = bsxfun(@times,kSpace,fftshift_mask); % fftshift on first 2 dimension can put the result image into right phase. Or the recon-image has to be multiplied by fftshift_mask to be in the right phase. Not doing so makes the image sharper, and the k-space has center in the center.
%phase_mod = double(logical(abs(sum(kSpace(:,:,:,1,1,1,:),7))));
%phase_mod(:,:,:,1,2) = exp(-1i*2*pi/3)*phase_mod(:,:,:,1,1);
%phase_mod(:,:,:,1,3) = exp(-1i*4*pi/3)*phase_mod(:,:,:,1,1);
%phase_mod_all = repmat(phase_mod,[sx,1,1,1,1]); clear phase_mod
%phase_mod = single(logical(abs(kSpace(:,:,:,1,1,1))));

%if nSMS>1
%    phase_mod(:,:,:,:,2:nSMS) = pre_interp_radial_to_cart(phase_mod_all(:,:,:,:,2:nSMS),x_coor,y_coor);
%end
    
%kSpace = fftshift(fftshift(kSpace,1),2);
%phase_mod = fftshift(fftshift(phase_mod,1),2);
para.CPUtime.prepare_kSpace = toc(t1);toc(t1);fprintf('\n');
  