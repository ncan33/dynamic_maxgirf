function [kSpace,phase_mod,N,para] = prepare_kSpace_NUFFT_ye(kSpace_all,theta_all,phase_mod_all,RayPosition,para)

[sx,~,nc,~,ns] = size(kSpace_all);
%ifasym              = para.asymmetry;
nSMS                = para.Recon.nSMS;
%ifNUFFT             = para.ifNUFFT;
kCenter             = para.kSpace_center;
%ifA3J               = para.FirstEstimationA3J;
selected_rays_start = para.selected_rays_start;
selected_rays_end   = para.selected_rays_end;

nor = selected_rays_end - selected_rays_start + 1; para.Recon.nor = nor;

if para.cSMS
    temp = find(RayPosition==1);
    RayPosition(temp(end):end) = 0; % delete the last few rays
    RayPosition(temp(1):temp(2)-1) = 0; %deleta first few rays
    
    index_1 = find(RayPosition==1);
    nor_all = diff(index_1);
    nof_low_ray_number = find(nor_all<selected_rays_end+1);
    
    for i=1:length(nof_low_ray_number)
        RayPosition(index_1(nof_low_ray_number(i)):index_1(nof_low_ray_number(i)+1)-1) = 0;
    end
    
    RayPosition(index_1(end):end) = 0;% delete last few rays again
end

RayIndex = logical((RayPosition > selected_rays_start-1) .* (RayPosition < selected_rays_end+1));
theta = theta_all(RayIndex);
nof = length(theta)/nor; para.Recon.nof = nof;
theta = reshape(theta,[1,nor,nof]);
kSpace = kSpace_all(:,RayIndex,:,:,:);
kSpace = reshape(kSpace,[sx,nor,nof,nc,1,ns]);
phase_mod = phase_mod_all(1,RayIndex,1,:);
phase_mod = reshape(phase_mod,[1,nor,nof,1,nSMS]);


disp('Creating NUFFT object...');tic   
[kx, ky] = get_k_coor(sx,theta,0,kCenter);
%coor = x_coor + 1i*y_coor;
    
%if ifA3J
%    [x_coor, y_coor] = get_k_coor(sx,theta_A3J,ifNUFFT,kCenter);
%    coor = x_coor + 1i*y_coor;
%end

%%%%% density compensation for NUFFT
%center_shift = kCenter - sx/2 -1;
%W = designFilter(sx,center_shift,'ram-lak');
%%%%% NUFFT part

%if ifasym
%    W(asymRayStart:asymRayEnd) = [];
%    coor(asymRayStart:asymRayEnd,:,:) = [];
%    kSpace(asymRayStart:asymRayEnd,:,:,:,:) = [];
%    if ifA3J
%        coor_FE(asymRayStart:asymRayEnd,:,:) = [];
%        kSpace_FE(asymRayStart:asymRayEnd,:,:,:,:) = [];
%    end
%end

switch para.Recon.interp_method
    case 'NUFFT'
        %[kx,ky] = GROG.Trajectory_correction_new(kSpace,kx,ky,500);
        N = NUFFT.init_new(kx,ky,para.over_sampling,para.core_size); %create the NUFFT structure
    case 'Toeplitz'
        if nSMS~=1
        
        N = cell(1,nSMS);
        phase_mod = round(imag(phase_mod(:,:,:,:,2)));
        SMS(1,:,:,1,1,1) = phase_mod==0;
        SMS(1,:,:,1,1,2) = phase_mod==1;
        SMS(1,:,:,1,1,3) = phase_mod==-1;
        sx_over = sx*para.over_sampling;
        
        kSpace_cart = single(zeros(sx_over,sx_over,nof,nc,ns,1,nSMS));
        mask_Toep = single(zeros(sx_over,sx_over,nof,1,ns,1,nSMS));
        for j=1:nSMS
            kx_temp = kx(repmat(SMS(:,:,:,:,:,j),[sx 1 1]));
            ky_temp = ky(repmat(SMS(:,:,:,:,:,j),[sx 1 1]));
            kx_temp = reshape(kx_temp,[sx nor/nSMS nof]);
            ky_temp = reshape(ky_temp,[sx nor/nSMS nof]);
            N{j} = NUFFT.init_new(kx_temp,ky_temp,para.over_sampling,para.core_size);
            
            kSpace_radial_temp = kSpace(repmat(SMS(:,:,:,:,:,j),[sx 1 1 nc]));
            kSpace_radial_temp = reshape(kSpace_radial_temp,[sx nor/nSMS nof nc]);
            kSpace_cart(:,:,:,:,:,:,j) = NUFFT.rad2cart(kSpace_radial_temp.*N{j}.W,N{j});
            mask_Toep(:,:,:,:,:,:,j) = NUFFT.rad2cart(bsxfun(@times,ones(size(kSpace_radial_temp(:,:,:,1))),N{j}.W),N{j});
        end
        
        else
            N = NUFFT.init_new(kx,ky,para.over_sampling,para.core_size);
            kSpace_cart = NUFFT.rad2cart(kSpace.*N.W,N);
            mask_Toep = NUFFT.rad2cart(N.W,N);
        end
        
        kSpace = kSpace_cart;
        phase_mod = abs(mask_Toep);
end

    
%if ifA3J
%    NUFFT_FE = FftTools.MultiNufft(coor_FE,W,0,[sx sx]); %create the NUFFT structure
%else
%    NUFFT_FE = NUFFT;
%    phase_mod_FE_conj = phase_mod_conj;
%    phase_mod_FE = phase_mod;
%    kSpace_FE = kSpace;
%end

%img_FE = NUFFT_FE'* bsxfun(@times,kSpace_FE,phase_mod_FE_conj); % iNUFFT
    
%if noSensMap
%    kSpace_SM = bsxfun(@times,kSpace_all,phase_mod_all);
%    [x_coor_SM, y_coor_SM] = get_k_coor(sx,theta_all,ifNUFFT,kCenter);
%    coor_SM = x_coor_SM + 1i*y_coor_SM;
%    if ifasym
%        coor_SM(asymRayStart:asymRayEnd,:,:) = [];
%        kSpace_SM(asymRayStart:asymRayEnd,:,:,:,:) = [];
%    end
%    NUFFT_SM = FftTools.MultiNufft(coor_SM,W,0,[sx sx]);
        %NUFFT_SM.kSize_z = no_comp;
%end
para.CPUtime.create_NUFFT_obj = toc;toc;
