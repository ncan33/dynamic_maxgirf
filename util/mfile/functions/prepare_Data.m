function [Data,para] = prepare_Data(kSpace_all,RayPosition,para)
% [kSpace,phase_mod,para] = prepare_kSpace(kSpace_all,theta_all,phase_mod_all,RayPosition,para)
%kSpace_all(1:152,:,:) = 0;
%kSpace_all(297:end,:,:) = 0;
t1 = tic;
[sx,~,no_comp,~,ns] = size(kSpace_all);
%nof                 = para.Recon.nof;
nSMS                = para.Recon.nSMS;
selected_rays_start = para.selected_rays_start;
selected_rays_end   = para.selected_rays_end;
%ifNUFFT             = para.ifNUFFT;
kCenter             = para.kSpace_center;
interp_method       = para.Recon.interp_method;

nor = selected_rays_end - selected_rays_start + 1; para.Recon.nor = nor;

if para.cSMS
    temp = find(RayPosition==1);
    RayPosition(temp(end-2):end) = 0; % delete the last few rays
    RayPosition(temp(1):temp(2)-1) = 0; %deleta first few rays
    
    index_1 = find(RayPosition==1);
    nor_all = diff(index_1);
    nof_low_ray_number = find(nor_all<selected_rays_end);
    
    for i=1:length(nof_low_ray_number)
        RayPosition(index_1(nof_low_ray_number(i)):index_1(nof_low_ray_number(i)+1)-1) = 0;
    end
end

RayIndex = logical((RayPosition > selected_rays_start-1) .* (RayPosition < selected_rays_end+1));
theta_all     = get_angle_mod(para);% radial sampling angle
phase_mod_all = single(get_phase_mod(para));% if SMS

theta = theta_all(RayIndex); nof = length(theta)/nor;
theta = reshape(theta,[1,nor,nof]);
kSpace_radial = kSpace_all(:,RayIndex,:,:,:);
kSpace_radial = reshape(kSpace_radial,[sx,nor,nof,no_comp,1,ns]);
phase_mod = phase_mod_all(1,RayIndex,1,:);
phase_mod = reshape(phase_mod,[1,nor,nof,1,nSMS]);

if ~isfield(para.Recon,'RF_frames') && isfield(para.Recon,'PD_frames')
    para.Recon.RF_frames = para.Recon.PD_frames(end)+1:nof;
end

%scale_kSpace = max(abs(kSpace_radial(:))); % scale to the max of kspace
%kSpace_radial = kSpace_radial/scale_kSpace*prod(para.Recon.image_size)*100; % try to scale image into order of~10
% be careful that single precision can only handle ~10e-7. if smaller than that it would be totally wrong, so careful about the scale! 

%%%%% pre-interpolation or NUFFT
disp('Pre-interpolate onto Cartesian space...')

matObj = matfile(para.dir.PCA_dir);% detect if there's trajectory
varlist = who(matObj,'Kx');
if para.asymmetry
    sx = para.AsymmetryRayEnd;
end
if ~isempty(varlist)
    load(para.dir.PCA_dir,'Kx','Ky')
    kx = double(squeeze(Kx(:,RayIndex))); ky = double(squeeze(Ky(:,RayIndex)));
    kx = reshape(kx,[sx,nor,nof]); ky = reshape(ky,[sx,nor,nof]);
else
    [kx, ky] = get_k_coor(sx,theta,0,kCenter);
end

varlist = who(matObj,'kSpace_info');
if ~isempty(varlist)
    load(para.dir.PCA_dir,'kSpace_info')
    if isfield(kSpace_info,'kx')
        kx = double(squeeze(kSpace_info.kx(:,RayIndex))); ky = double(squeeze(kSpace_info.ky(:,RayIndex)));
        kx = reshape(kx,[sx,nor,nof]); ky = reshape(ky,[sx,nor,nof]);
%         theta = atan(ky(1,:,:)./kx(1,:,:));
%         kx = kx - cos(theta)*0.5;
%         ky = ky - sin(theta)*0.5;
    end
else
    [kx, ky] = get_k_coor(sx,theta,0,kCenter);
end

switch interp_method
    case 'grid3'
        if nSMS ~= 1
            phase_mod = round(imag(phase_mod(:,:,:,:,2)));
            SMS(1,:,:,1,1,1) = phase_mod==0;
            SMS(1,:,:,1,1,2) = phase_mod==1;
            SMS(1,:,:,1,1,3) = phase_mod==-1;
            sx_over = sx*para.over_sampling;
            Data.kSpace = single(zeros(sx_over,sx_over,nof,no_comp,ns,1,nSMS));
            for j=1:nSMS
                kx_temp = kx(repmat(SMS(:,:,:,:,:,j),[sx 1 1]));
                ky_temp = ky(repmat(SMS(:,:,:,:,:,j),[sx 1 1]));
                kx_temp = reshape(kx_temp,[sx nor/nSMS nof]);
                ky_temp = reshape(ky_temp,[sx nor/nSMS nof]);
                kSpace_radial_temp = kSpace_radial(repmat(SMS(:,:,:,:,:,j),[sx 1 1 no_comp]));
                kSpace_radial_temp = reshape(kSpace_radial_temp,[sx nor/nSMS nof no_comp]);
                Data.kSpace(:,:,:,:,:,:,j) = pre_interp_radial_to_cart_over(kSpace_radial_temp,kx_temp,ky_temp,para.over_sampling);
            end
            para.Recon.type = 'seperate SMS';
        else
            Data.kSpace = pre_interp_radial_to_cart_oversample(kSpace_radial,kx,ky,para.over_sampling);
            para.Recon.type = '2D';
        end
        
    case'nn'
        N = NN.NN_init(kx,ky,para.over_sampling);
        Data.kSpace = NN.NN_rad2cart(kSpace_radial,N);
        para.Recon.type = '2D';

    case {'GROG' 'GNUFFT'}

        [Data.G,Data.kSpace] = GROG.GROG_seperate_SMS_GNUFFT(squeeze(kSpace_radial),kx,ky,phase_mod, para);

        if nSMS == 1            
            para.Recon.type = '2D';
        else
            para.Recon.type = 'seperate SMS';
        end
        
        %s = round(Data.G{1}.core_size/2);
        %Data.kSpace = Data.kSpace(s+1:s+Data.G{1}.sx_over, s+1:s+Data.G{1}.sx_over,:,:,:,:,:);

        switch para.Recon.interp_method
            case 'GNUFFT'
                para.Recon.type = 'GNUFFT';
        end

    case 'NUFFT'
        Data.kSpace = kSpace_radial;
        Data.N = NUFFT.init_new(kx,ky,para.over_sampling,para.core_size);
        im = NUFFT.NUFFT_adj_new(kSpace_radial.*conj(phase_mod),Data.N);
        if nSMS == 1
            Data.sens_map = get_sens_map(im,'2D');
        else
            Data.sens_map = get_sens_map(im,'SMS');
            Data.phase_mod = phase_mod;
        end
        Data.first_est = bsxfun(@times,im,conj(Data.sens_map));
        Data.first_est = sum(Data.first_est,4);
        %para.Recon.type = 'NUFFT';
        %if size(Data.first_est,3)>100
            para.Recon.type = 'NUFFT';
        %end
        scale_image = max(abs(Data.first_est(:)));
        para.Recon.weight_tTV = scale_image*para.weight_tTV;
        para.Recon.weight_sTV = scale_image*para.weight_sTV;

        return
        
    case 'Toeplitz'
        if nSMS~=1
            N = cell(1,nSMS);
            phase_mod = round(imag(phase_mod(:,:,:,:,2)));
            SMS(1,:,:,1,1,1) = phase_mod==0;
            SMS(1,:,:,1,1,2) = phase_mod==1;
            SMS(1,:,:,1,1,3) = phase_mod==-1;
            sx_over = sx*para.over_sampling;
            para.Recon.kSpace_size = [sx_over,sx_over];
            
            Data.kSpace = single(zeros(sx_over,sx_over,nof,no_comp,ns,1,nSMS));
            Data.mask = single(zeros(sx_over,sx_over,nof,1,ns,1,nSMS));
            for j=1:nSMS
                kx_temp = kx(repmat(SMS(:,:,:,:,:,j),[sx 1 1]));
                ky_temp = ky(repmat(SMS(:,:,:,:,:,j),[sx 1 1]));
                kx_temp = reshape(kx_temp,[sx nor/nSMS nof]);
                ky_temp = reshape(ky_temp,[sx nor/nSMS nof]);
                N{j} = NUFFT.init_new(kx_temp,ky_temp,para.over_sampling,para.core_size);
                
                kSpace_radial_temp = kSpace_radial(repmat(SMS(:,:,:,:,:,j),[sx 1 1 no_comp]));
                kSpace_radial_temp = reshape(kSpace_radial_temp,[sx nor/nSMS nof no_comp]);
                Data.kSpace(:,:,:,:,:,:,j) = NUFFT.rad2cart(kSpace_radial_temp.*N{j}.W,N{j});
                Data.mask(:,:,:,:,:,:,j) = NUFFT.rad2cart(bsxfun(@times,ones(size(kSpace_radial_temp(:,:,:,1))),N{j}.W),N{j});
            end
            
            para.Recon.type = 'Toeplitz SMS';
            %{
            if nof>50
                para.Recon.type = 'seperate SMS less memory';
            else
                para.Recon.type = 'seperate SMS';
            end
            %}
            Data.Apodizer = N{1}.Apodizer;
            
            kSpace_sms(:,:,:,:,1) = sum(Data.kSpace,7);
            kSpace_sms(:,:,:,:,2) = Data.kSpace(:,:,:,:,:,:,1) + exp(-1i*2*pi/3)*Data.kSpace(:,:,:,:,:,:,2) + exp(-1i*4*pi/3)*Data.kSpace(:,:,:,:,:,:,3);
            kSpace_sms(:,:,:,:,3) = Data.kSpace(:,:,:,:,:,:,1) + exp(-1i*4*pi/3)*Data.kSpace(:,:,:,:,:,:,2) + exp(-1i*2*pi/3)*Data.kSpace(:,:,:,:,:,:,3);
            
            im = ifft2(kSpace_sms);
            im(para.Recon.sx+1:end,:,:,:,:,:) = [];
            im(:,para.Recon.sx+1:end,:,:,:,:) = [];
            
            im = bsxfun(@times,im,N{1}.Apodizer);
            sens_map = get_sens_map(im,'SMS');
        else
            N = NUFFT.init_new(kx,ky,para.over_sampling,para.core_size);
            Data.kSpace = NUFFT.rad2cart(kSpace_radial.*N.W,N);
            Data.mask = NUFFT.rad2cart(N.W,N);
            Data.Apodizer = N.Apodizer;
            para.Recon.type = '2D';
            
            im = ifft2(Data.kSpace);
            im(para.Recon.sx+1:end,:,:,:,:,:) = [];
            im(:,para.Recon.sx+1:end,:,:,:,:) = [];
            im = bsxfun(@times,im,N.Apodizer);
            sx_over = size(Data.kSpace,1);
            para.Recon.kSpace_size = [sx_over,sx_over];
            sens_map = get_sens_map(im,'2D');
            
        end
        Data.mask = abs(Data.mask);
        Data.sens_map = sens_map;
        Data.first_est = single(sum(bsxfun(@times, im, conj(sens_map)),4));
        
        if para.image_orintation == 0
            para.image_orintation = orintation_detection(abs(fftshift(ifft2(sum(sum(sum(Data.kSpace,3),4),7)))));
            Data.kSpace = orintate_image(Data.kSpace,para.image_orintation);
            Data.mask = orintate_image(Data.mask,para.image_orintation);
            Data.sens_map = orintate_image(Data.sens_map,para.image_orintation);
            Data.first_est = orintate_image(Data.first_est,para.image_orintation);
        else
            Data.kSpace = orintate_image(Data.kSpace,para.image_orintation);
            Data.mask = orintate_image(Data.mask,para.image_orintation);
            Data.sens_map = orintate_image(Data.sens_map,para.image_orintation);
            Data.first_est = orintate_image(Data.first_est,para.image_orintation);
        end
        
        return
        
end

% zero pad
if isfield(para.Recon,'zero_pad')
    zero_pad = para.Recon.zero_pad;
    para.Recon.sx = zero_pad;
    para.Recon.image_size(1) = zero_pad;
    para.Recon.image_size(2) = zero_pad;
    
    kSpace = zeros(para.Recon.sx,para.Recon.sx,para.Recon.nof,para.Recon.no_comp,1,1,3,'single');
    center = (para.Recon.sx-sx)/2+1:(para.Recon.sx-sx)/2+sx;
    kSpace(center,center,:,:,:,:,:) = Data.kSpace;
    Data.kSpace = kSpace;
end

%low res
if isfield(para.Recon,'low_res')
    Data.kSpace = crop_half_FOV(Data.kSpace);
    para.Recon.image_size(1) = size(Data.kSpace,1);
    para.Recon.image_size(2) = size(Data.kSpace,2);
    
end
% orintation
if para.image_orintation == 0
    para.image_orintation = orintation_detection(abs(fftshift(ifft2(sum(sum(sum(Data.kSpace,3),4),7)))));
    Data.kSpace = orintate_image(Data.kSpace,para.image_orintation);
else
    Data.kSpace = orintate_image(Data.kSpace,para.image_orintation);
end

Data.kSpace = fftshift2(Data.kSpace);
Data.mask = logical(abs(Data.kSpace(:,:,:,1,:,:,:)));

switch para.Recon.type
    case 'seperate SMS'
        kSpace_sms(:,:,:,:,1) = sum(Data.kSpace,7);
        kSpace_sms(:,:,:,:,2) = Data.kSpace(:,:,:,:,:,:,1) + exp(-1i*2*pi/3)*Data.kSpace(:,:,:,:,:,:,2) + exp(-1i*4*pi/3)*Data.kSpace(:,:,:,:,:,:,3);
        kSpace_sms(:,:,:,:,3) = Data.kSpace(:,:,:,:,:,:,1) + exp(-1i*4*pi/3)*Data.kSpace(:,:,:,:,:,:,2) + exp(-1i*2*pi/3)*Data.kSpace(:,:,:,:,:,:,3);

        Data.kSpace = ifft2(Data.kSpace);
        Data.kSpace = fftshift2(Data.kSpace);
        Data.kSpace = fft2(Data.kSpace);
        Data.kSpace = Data.kSpace .* Data.mask;
        %Data.kSpace = fft2(fftshift2(ifft2(Data.kSpace))).*Data.mask;
        im = ifft2(kSpace_sms);
        im = fftshift2(im);
        kSpace_sms = fft2(im).*logical(sum(Data.mask,7));
        sx_over = size(im,1);
        cut = (sx_over - para.Recon.image_size(1))/2;
        para.Recon.kSpace_size = [sx_over,sx_over];
        im([1:cut,end-cut+1:end],:,:,:,:) = [];
        im(:,[1:cut,end-cut+1:end],:,:,:) = [];
        Data.sens_map = get_sens_map(im,'SMS');
        Data.filter = ramp_filter_for_pre_interp(para);
        
        im = ifft2(kSpace_sms.*Data.filter);
        im([1:cut,end-cut+1:end],:,:,:,:) = [];
        im(:,[1:cut,end-cut+1:end],:,:,:) = [];
        Data.first_est = single(sum(bsxfun(@times, im, conj(Data.sens_map)),4));
%         if nof>1
%             para.Recon.type = 'seperate SMS less memory';
%         end
    case '2D'
        im = ifft2(Data.kSpace);
        im = fftshift2(im);
        Data.kSpace = fft2(im).*Data.mask;
        sx_over = size(im,1);
        cut = (sx_over - para.Recon.image_size(1))/2;
        para.Recon.kSpace_size = [sx_over,sx_over];
        im([1:cut,end-cut+1:end],:,:,:,:) = [];
        im(:,[1:cut,end-cut+1:end],:,:,:) = [];
        Data.sens_map = get_sens_map(im,'2D');
        Data.filter = ramp_filter_for_pre_interp(para);
        
        im = ifft2(Data.kSpace.*Data.filter);
        im([1:cut,end-cut+1:end],:,:,:,:) = [];
        im(:,[1:cut,end-cut+1:end],:,:,:) = [];
        Data.first_est = single(sum(bsxfun(@times, im, conj(Data.sens_map)),4));
end

scale_image = max(abs(Data.first_est(:)));
para.Recon.weight_tTV = scale_image*para.weight_tTV;
para.Recon.weight_sTV = scale_image*para.weight_sTV;

para.CPUtime.prepare_kSpace = toc(t1);toc(t1);fprintf('\n');

