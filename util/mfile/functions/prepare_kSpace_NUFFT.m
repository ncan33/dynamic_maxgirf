function [kSpace,phase_mod,NUFFT,para] = prepare_kSpace_NUFFT(kSpace_all,theta_all,phase_mod_all,RayPosition,para)

[sx,~,no_comp,~,ns] = size(kSpace_all);
ifasym              = para.asymmetry;

nSMS                = para.Recon.nSMS;
ifNUFFT             = para.ifNUFFT;
kCenter             = para.kSpace_center;
%ifA3J               = para.FirstEstimationA3J;
selected_rays_start = para.selected_rays_start;
selected_rays_end   = para.selected_rays_end;

nor = selected_rays_end - selected_rays_start + 1; para.Recon.nor = nor;
RayIndex = logical((RayPosition > selected_rays_start-1) .* (RayPosition < selected_rays_end+1));
theta = theta_all(RayIndex);
nof = length(theta)/nor; para.Recon.nof = nof;
theta = reshape(theta,[1,nor,nof]);
kSpace = kSpace_all(:,RayIndex,:,:,:);
kSpace = reshape(kSpace,[sx,nor,nof,no_comp,1,ns]);
phase_mod = phase_mod_all(1,RayIndex,1,:);
phase_mod = reshape(phase_mod,[1,nor,nof,1,nSMS]);


disp('Creating NUFFT object...');tic   
[x_coor, y_coor] = get_k_coor(sx,theta,ifNUFFT,kCenter);
coor = x_coor + 1i*y_coor;
    
%if ifA3J
%    [x_coor, y_coor] = get_k_coor(sx,theta_A3J,ifNUFFT,kCenter);
%    coor = x_coor + 1i*y_coor;
%end

%%%%% density compensation for NUFFT
center_shift = kCenter - sx/2 -1;
W = designFilter(sx,center_shift,'ram-lak');
%%%%% NUFFT part

if ifasym
    W(asymRayStart:asymRayEnd) = [];
    coor(asymRayStart:asymRayEnd,:,:) = [];
    kSpace(asymRayStart:asymRayEnd,:,:,:,:) = [];
%    if ifA3J
%        coor_FE(asymRayStart:asymRayEnd,:,:) = [];
%        kSpace_FE(asymRayStart:asymRayEnd,:,:,:,:) = [];
%    end
end

NUFFT = FftTools.MultiNufft(coor,W,0,[sx sx]); %create the NUFFT structure
    
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
