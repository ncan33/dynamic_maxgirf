function [Data,para] = prepare_Data_3D_new(kSpace_radial,kSpace_info,para)
disp('Pre-interpolate into Cartesian space...');t1=tic;
if length(size(kSpace_radial))<5
    kSpace_radial = permute(kSpace_radial,[1,2,4,5,3]);
end
slice_idx = kSpace_info.slice_idx;
if min(slice_idx(:)) <= 0
    slice_idx = slice_idx - min(slice_idx) + 1;
end
theta = kSpace_info.angle_mod;
frames = kSpace_info.frames;

nos = max(slice_idx(:)) + 2;
%nof = max(frames(:));
sx = size(kSpace_radial,1);
nc = size(kSpace_radial,5);
nor = diff(find([1,diff(frames),length(slice_idx)]));
nof = length(nor);
[kx_all,ky_all] = get_k_coor(sx,theta,0,round(sx/2)+1);
kSpace = zeros(sx,1,nos,nof,nc,'single');
kx = zeros(sx,1,nos,nof,'single');
ky = zeros(sx,1,nos,nof,'single');
for i=1:nof
    rays = (1:nor(i)) + sum(nor(1:i-1));
    kSpace_temp = squeeze(kSpace_radial(:,rays,:,:,:));
    slice_idx_temp = slice_idx(rays);
    kx_temp = kx_all(:,rays);
    ky_temp = ky_all(:,rays);
    for j=1:nos-2
        rays_slice = slice_idx_temp==j;
        kSpace(:,1:sum(rays_slice),j,i,:) = kSpace_temp(:,rays_slice,:);
        kx(:,1:sum(rays_slice),j,i) = kx_temp(:,rays_slice);
        ky(:,1:sum(rays_slice),j,i) = ky_temp(:,rays_slice);
    end
end
kSpace_radial = kSpace;

nor = kSpace_radial(para.kSpace_center,:,:,1,1);
nor = logical(abs(nor));
para.Recon.nor = squeeze(sum(nor,2));
para.Recon.sx = size(kSpace_radial,1)*para.over_sampling;

%scale_kSpace = max(abs(kSpace_radial(:))); % scale to the max of kspace
%kSpace_radial = kSpace_radial/scale_kSpace*para.Recon.sx^2*800; % try to scale image into order of~10

switch para.Recon.interp_method
    case 'GROG'
        Data.kSpace = GROG.GROG_3D(kSpace_radial,kx,ky,0,1);
    case 'Toeplitz'
        [Data.kSpace,Data.mask,Data.Apodizer] = Toeplitz_3D(kSpace_radial,kx,ky,para);
    case 'NUFFT'
        Data = NUFFT.ThreeD_init(kSpace_radial,kx,ky,para);
    case 'grid3'
        [sx,nor,sz,nof,no_comp] = size(kSpace_radial);
        data = double(reshape(kSpace_radial,[sx*nor,sz*nof*no_comp]));
        kx = double(kx);
        ky = double(ky);

        x = reshape(kx,[sx*nor,sz*nof]);
        y = reshape(ky,[sx*nor,sz*nof]);

        x = repmat(x,[1 no_comp]);
        y = repmat(y,[1 no_comp]);
        
        Xr = round(x);
        Yr = round(y);

        kSpace_cart = single(zeros((sx+1)*(sx+1),sz*nof*no_comp));
        kSpace_r = single(zeros(sx*nor,sz*nof*no_comp));
        
        for i=1:sz*nof*no_comp
            warning off
            index = data(:,i) ~= 0;
            kSpace_r(index,i) = griddata(x(index,i),y(index,i),data(index,i),Xr(index,i),Yr(index,i));
        end
        
        kSpace_r(isnan(kSpace_r)) = 0;
        
        indx = sub2ind([sx+1,sx+1,sz*nof*no_comp],Xr+sx/2+1,Yr+sx/2+1);

        for i=1:sz*nof*no_comp
            kSpace_cart(indx(:,i),i) = kSpace_r(:,i);
        end
        
        kSpace_cart = reshape(kSpace_cart,[sx+1,sx+1,sz,nof,no_comp]);
        kSpace_cart(1,:,:,:,:,:) = [];
        kSpace_cart(:,1,:,:,:,:) = [];
        Data.kSpace = kSpace_cart;
end
%save(interp_dir,'kSpace_cart','-v7.3')

switch para.Recon.interp_method
    case 'NUFFT'
        im = NUFFT.NUFFT_adj_new(Data.kSpace,Data.N);
        Data.sens_map = get_sens_map(im,'3D');
        para.Recon.type = '3D NUFFT';
    otherwise
        if para.image_orintation == 0
            para.image_orintation = orintation_detection(abs(fftshift(ifft2(sum(sum(sum(Data.kSpace,3),4),7)))));
            Data.kSpace = orintate_image(Data.kSpace,para.image_orintation);
        else
            Data.kSpace = orintate_image(Data.kSpace,para.image_orintation);
        end

        Data.kSpace(isnan(Data.kSpace)) = 0;
        Data.kSpace = fftshift3(Data.kSpace);
        
        Data.mask = logical(abs(Data.kSpace(:,:,:,:,1)));
        
        Data.kSpace = ifft3(Data.kSpace);
        Data.kSpace = fftshift3(Data.kSpace);
        % move the first slice to end to ordering
        %Data.kSpace = circshift(Data.kSpace,-1,3);

        Data.kSpace = fft3(Data.kSpace);
        Data.kSpace = Data.kSpace.*Data.mask;
        im = ifft3(Data.kSpace);

        Data.filter = ramp_filter_for_pre_interp_3D(para);
        
        Data.sens_map = get_sens_map(im,'3D');
        Data.first_est = sum(bsxfun(@times,ifft3(Data.kSpace.*Data.filter),conj(Data.sens_map)),5);
        para.Recon.type = '3D';
end
%%%%% sensitivity map and first estimation
%[sens_map,para] = get_sens_map_3D(kSpace,para);
%Data.kSpace = Data.kSpace(:,:,:,para.SelectedFrames,:);
%Data.mask = Data.mask(:,:,:,para.SelectedFrames,:);
%Data.first_est = Data.first_est(:,:,:,para.SelectedFrames,:);

para.Recon.no_comp = size(Data.kSpace,5);
para.CPUtime.prepare_kSpace = toc(t1);toc(t1);fprintf('\n');
end