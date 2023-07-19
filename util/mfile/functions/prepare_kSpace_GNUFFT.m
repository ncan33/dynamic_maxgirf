function [G,kSpace,kSpace_radial,para] = prepare_kSpace_GNUFFT(kSpace_all,theta_all,phase_mod_all,RayPosition,para)
% [kSpace,phase_mod,para] = prepare_kSpace(kSpace_all,theta_all,phase_mod_all,RayPosition,para)

t1 = tic;
[sx,~,no_comp,~,ns] = size(kSpace_all);
%nof                 = para.Recon.nof;
nSMS                = para.Recon.nSMS;
selected_rays_start = para.selected_rays_start;
selected_rays_end   = para.selected_rays_end;
ifNUFFT             = para.ifNUFFT;
kCenter             = para.kSpace_center;
interp_method       = para.Recon.interp_method;

nor = selected_rays_end - selected_rays_start + 1; para.Recon.nor = nor;

if para.cSMS
    temp = find(RayPosition==1);
    RayPosition(temp(end-2):end) = 0; % delete the last few rays
    RayPosition(temp(1):temp(2)-1) = 0; %deleta first few rays
    
    index_1 = find(RayPosition==1);
    nor_all = diff(index_1);
    nof_low_ray_number = find(nor_all<220);
    
    for i=1:length(nof_low_ray_number)
        RayPosition(index_1(nof_low_ray_number(i)):index_1(nof_low_ray_number(i)+1)-1) = 0;
    end
end

RayIndex = logical((RayPosition > selected_rays_start-1) .* (RayPosition < selected_rays_end+1));
theta = theta_all(RayIndex); nof = length(theta)/nor;
theta = reshape(theta,[1,nor,nof]);
kSpace = kSpace_all(:,RayIndex,:,:,:);
kSpace = reshape(kSpace,[sx,nor,nof,no_comp,1,ns]);
phase_mod = phase_mod_all(1,RayIndex,1,:);
phase_mod = reshape(phase_mod,[1,nor,nof,1,nSMS]);

%%%%% pre-interpolation
disp('Pre-interpolate into Cartesian space...')

matObj = matfile(para.dir.PCA_dir);
varlist = who(matObj,'Kx');
if para.asymmetry
    sx = para.AsymmetryRayEnd;
end
if ~isempty(varlist)
    load(para.dir.PCA_dir,'Kx','Ky')
    x_coor = double(squeeze(Kx)); y_coor = double(squeeze(Ky));
else
    [x_coor, y_coor] = get_k_coor(sx,theta,ifNUFFT,kCenter);
end

switch interp_method
    case 'grid3'
        %kSpace = pre_interp_radial_to_cart(kSpace,x_coor,y_coor);
        if nSMS ~= 1
            phase_mod = round(imag(phase_mod(:,:,:,:,2)));
            SMS(1,:,:,1,1,1) = phase_mod==0;
            SMS(1,:,:,1,1,2) = phase_mod==1;
            SMS(1,:,:,1,1,3) = phase_mod==-1;
            sx_over = sx*para.over_sampling;
            kSpace_cart = single(zeros(sx_over,sx_over,nof,no_comp,ns,1,nSMS));
            for j=1:nSMS
                kx_temp = x_coor(repmat(SMS(:,:,:,:,:,j),[sx 1 1]));
                ky_temp = y_coor(repmat(SMS(:,:,:,:,:,j),[sx 1 1]));
                kx_temp = reshape(kx_temp,[sx nor/nSMS nof]);
                ky_temp = reshape(ky_temp,[sx nor/nSMS nof]);
                
                kSpace_radial_temp = kSpace(repmat(SMS(:,:,:,:,:,j),[sx 1 1 no_comp]));
                kSpace_radial_temp = reshape(kSpace_radial_temp,[sx nor/nSMS nof no_comp]);
                kSpace_cart(:,:,:,:,:,:,j) = pre_interp_radial_to_cart_over(kSpace_radial_temp,kx_temp,ky_temp,para.over_sampling);
            end
            kSpace = kSpace_cart;
        else
            kSpace = pre_interp_radial_to_cart_oversample(kSpace,x_coor,y_coor,para.over_sampling);
        end
        G = 0;
        kSpace_radial = 0;
        
    case'nn'
        
        G = NN.NN_init(kSpace,x_coor,y_coor);
        kSpace_radial = kSpace;
        kSpace = NN.NN_rad2cart(kSpace,G);
        %kSpace = permute(kSpace,[1 2 4 3 5 6]);
        %kSpace_c = zeros([sx sx no_comp nof 1 ns]);
        %parfor i=1:ns
        %    kSpace_c(:,:,:,:,1,i) = interp_nn(kSpace(:,:,:,:,1,i),x_coor, y_coor);
        %end
        %kSpace = permute(kSpace_c,[1 2 4 3 5 6]); clear kSpace_c
    case {'GROG' 'GNUFFT'}
        if nSMS == 1
            [G,kSpace,kSpace_radial] = GROG.GROG_seperate_SMS_GNUFFT(squeeze(kSpace),x_coor,y_coor,phase_mod, 0, para.np, para);
        else
            [G,kSpace,kSpace_radial] = GROG.GROG_seperate_SMS_GNUFFT(squeeze(kSpace),x_coor,y_coor,phase_mod(:,:,:,:,2), 0, para.np, para);
        end
        
        s = round(G{1}.core_size/2);
        kSpace = kSpace(s+1:s+G{1}.sx_over, s+1:s+G{1}.sx_over,:,:,:,:,:);
        %kSpace(1,:,:,:,:,:,:) = [];
        %kSpace(end,:,:,:,:,:,:) = [];
        %kSpace(:,1,:,:,:,:,:) = [];
        %kSpace(:,end,:,:,:,:,:) = [];

%        [x_all, y_all] = get_k_coor(sx,theta_all,ifNUFFT,kCenter);
%        x_all = reshape(x_all,[288 340 70]); y_all = reshape(y_all,[288 340 70]);
%        kSpace = GROG.GROG_seperate_SMS_train_core_using_all_data(reshape(kSpace_all,[288 340 70 8]),squeeze(kSpace),x_coor,y_coor,x_all,y_all, 0, 1);

    case 'NUFFT'
        
end

%fftshift_mask = single(ones(sx+1));
%fftshift_mask(1:2:sx+1,:) = -fftshift_mask(1:2:sx+1,:);
%fftshift_mask(:,1:2:sx+1) = -fftshift_mask(:,1:2:sx+1);
%kSpace = bsxfun(@times,kSpace,fftshift_mask); % fftshift on first 2 dimension can put the result image into right phase. Or the recon-image has to be multiplied by fftshift_mask to be in the right phase. Not doing so makes the image sharper, and the k-space has center in the center.

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
  