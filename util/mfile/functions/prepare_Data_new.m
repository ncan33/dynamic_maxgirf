function [Data,para] = prepare_Data_new(kSpace_all,kSpace_info,para)
% [kSpace,phase_mod,para] = prepare_kSpace(kSpace_all,theta_all,phase_mod_all,RayPosition,para)
%kSpace_all(1:152,:,:) = 0;
%kSpace_all(297:end,:,:) = 0;        
t1 = tic;

selected_rays_start = para.selected_rays_start;
selected_rays_end   = para.selected_rays_end;
nSMS                = para.Recon.nSMS;
nr = selected_rays_end - selected_rays_start + 1;

theta = kSpace_info.angle_mod;
phase = kSpace_info.phase_mod;
frames = kSpace_info.frames;

nof = max(frames(:));
sx = size(kSpace_all,1);
nc = size(kSpace_all,5);
nor = diff(find([1,diff(frames),length(theta)]));
kSpace_radial = zeros(sx,1,nof,nc,'single');
theta_all = zeros(1,1,nof,'single');
phase_all = repmat([0,1,2],[1,nr/nSMS,nof]);
for i=1:nof
    rays = (1:nor(i)) + sum(nor(1:i-1));
    kSpace_temp = squeeze(kSpace_all(:,rays,:,:,:));
    theta_temp = theta(rays);
    phase_temp = phase(rays);
    rays = 1:nor(i);
    rays_idx = (rays >= selected_rays_start) & (rays <= selected_rays_end);
    nor_temp = floor(sum(rays_idx)/3);
    rays_idx(nor_temp*3+selected_rays_start:end) = false;
    kSpace_radial(:,1:sum(rays_idx),i,:) = kSpace_temp(:,rays_idx,:);
    theta_all(:,1:sum(rays_idx),i) = theta_temp(:,rays_idx,:);
    phase_all(:,1:sum(rays_idx),i) = phase_temp(:,rays_idx,:);
end

para.phase_mod = phase_all(:);
para.Recon.sy = nr;
para.Recon.sz = nof;
para.Recon.nor = nr;
para.Recon.sx = sx;
para.Recon.no_comp = nc;
phase_mod = get_phase_mod(para);
phase_mod = reshape(phase_mod,[1,nr,nof,1,nSMS]);

[sx,~,no_comp,~,ns] = size(kSpace_all);
%nof                 = para.Recon.nof;
kCenter             = para.kSpace_center;
interp_method       = para.Recon.interp_method;

%%%%% pre-interpolation or NUFFT
disp('Pre-interpolate onto Cartesian space...')

[kx, ky] = get_k_coor(sx,theta_all,0,kCenter);
para.Recon.image_size = [sx,sx];

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
        N = NN.NN_init(kSpace_radial,kx,ky);
        Data.kSpace = NN.NN_rad2cart(kSpace_radial,N);
        para.Recon.type = '2D';

    case {'GROG' 'GNUFFT'}

        [Data.G,Data.kSpace,~] = GROG.GROG_seperate_SMS_GNUFFT(squeeze(kSpace_radial),kx,ky,phase_mod, para);

        if nSMS == 1            
            para.Recon.type = '2D';
        else
            para.Recon.type = 'seperate SMS';
        end
        
        s = round(Data.G{1}.core_size/2);
%        Data.kSpace = Data.kSpace(s+1:s+Data.G{1}.sx_over, s+1:s+Data.G{1}.sx_over,:,:,:,:,:);
        
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
        
        para.CPUtime.prepare_kSpace = toc(t1);toc(t1);fprintf('\n');

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

% orintation image
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
        %if nof>1
            %para.Recon.type = 'seperate SMS less memory';
        %end
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

