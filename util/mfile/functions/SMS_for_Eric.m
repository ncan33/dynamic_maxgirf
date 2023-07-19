addpath /v/raid1b/ytian/mfile/
addpath /v/raid1b/ytian/mfile/functions/

clear
im_size        = 288; %image size
para.Recon.sy  = 1000;%number of radial rays
para.Recon.sz  = 1;
para.dataAngle = 180;

para.phase_mod = 1;
para.angle_mod = 2;

theta = get_angle_mod(para);
phase_mod = repmat([1 -1],[1,para.Recon.sy/2]);
[kx,ky] = get_k_coor(288,theta,0,144);
im1 = phantom(im_size);
im2 = rot90(im1);

N = NUFFT.init(kx,ky,1.5,[6 6]);

k1 = NUFFT.NUFFT(im1,N);
k2 = NUFFT.NUFFT(im2,N);
k = k1+k2.*phase_mod;
out1 = NUFFT.NUFFT_adj(k1,N);
out2 = NUFFT.NUFFT_adj(k1.*phase_mod,N);
out3 = NUFFT.NUFFT_adj(k,N);
out4 = NUFFT.NUFFT_adj(k.*-phase_mod,N);

figure,
subplot(1,4,1)
imagesc(abs(out1))
axis image
axis off
colorbar

subplot(1,4,2)
imagesc(abs(out2))
axis image
axis off
colorbar

subplot(1,4,3)
imagesc(abs(out3))
axis image
axis off
colorbar

subplot(1,4,4)
imagesc(abs(out4))
axis image
axis off
colorbar

