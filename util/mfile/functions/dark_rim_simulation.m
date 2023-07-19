
ccc
nor = 120;
nof = 10;

im = false(256,256);
im(128,128) = true;
im = bwdist(im);
im = (im<50)*10+(im<75)*5;
im = repmat(im,[1,1,nof]);

theta = golden_angle(1:nor*nof);
theta = reshape(theta,[1,nor,nof]);
[kx,ky] = get_k_coor(256,theta,0,129);
N = NUFFT.init_new(kx,ky,1.5,[6,6]);
kSpace = NUFFT.NUFFT_new(im,N);
im_recon = NUFFT.NUFFT_adj_new(kSpace,N);

N.S = gpuArray(N.S);
kSpcae = gpuArray(kSpace);
im_recon = gpuArray(im_recon);

for i=1:500
    update = NUFFT.NUFFT_adj_new(kSpace - NUFFT.NUFFT_new(im_recon,N),N);
    update = update * 0.8 + compute_tTV_yt(im_recon,0.1,eps('single'));
    update = update + compute_sTV_yt(im_recon,0.0005,eps('single'));
    im_recon = im_recon + update;
end


figure
imagesc(abs(im_recon(:,:,5)))
colormap gray
axis image

figure
plot(abs(im_recon(128,:,5)))
hold on
plot(im(128,:,5))