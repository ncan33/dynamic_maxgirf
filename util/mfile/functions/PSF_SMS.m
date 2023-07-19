clear

para.Recon.sy=24;
para.Recon.sz=1;
para.dataAngle=180;

para.phase_mod = 1;
para.angle_mod = 2;

theta = get_angle_mod(para);
phase_mod = get_phase_mod(para);
[kx,ky] = get_k_coor(288,theta,0,144);
im = zeros(288);
im(144,144) = 1;
N = NUFFT.init(kx,ky,1.5,[6 6]);

k = NUFFT.NUFFT(im,N);
PSF = NUFFT.NUFFT_adj(bsxfun(@times,k,phase_mod),N);

figure,
subplot(1,2,1)
imagesc(abs(PSF(:,:,1,1)))
axis image
axis off
colorbar
subplot(1,2,2)
imagesc(abs(PSF(:,:,1,2)))
axis image
axis off
colorbar

max(max(abs(PSF(:,:,1,2))))/max(max(abs(PSF(:,:,1,1))))
