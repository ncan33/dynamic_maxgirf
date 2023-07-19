function [kSpace_cart,Toeplitz_mask,Apodizer] = Toeplitz_3D(kSpace_radial,kx,ky,para)
over_sampling = para.over_sampling;
kernel_size = para.kernel_size;

[nx,ny,nz,nof,nc] = size(kSpace_radial);
kSpace_cart = zeros([nx*over_sampling,nx*over_sampling,nz,nof,nc],'single');
Toeplitz_mask = zeros([nx*over_sampling,nx*over_sampling,nz,nof],'single');
for i=1:nz-2
    kSpace_temp = squeeze(kSpace_radial(:,:,i,:,:));
    kx_temp = squeeze(kx(:,:,i,:));
    ky_temp = squeeze(ky(:,:,i,:));
    
    mask = kSpace_temp == 0;
    mask = mask(144,:,1,1);
    
    kSpace_temp(:,mask,:,:) = [];
    kx_temp(:,mask,:,:) = [];
    ky_temp(:,mask,:,:) = [];
    
    N = NUFFT.init_new(kx_temp,ky_temp,over_sampling,kernel_size);
    kSpace_cart(:,:,i,:,:) = NUFFT.rad2cart(bsxfun(@times,kSpace_temp,N.W),N);
    Toeplitz_mask(:,:,i,:) = NUFFT.rad2cart(N.W,N);
end
Apodizer = N.Apodizer;