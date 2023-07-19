clear, close all
%{
im = zeros(288,288);
[X,Y] = meshgrid(1:288,1:288);
X = X - 144.5;
Y = Y - 144.5;
heart = sqrt(X.^2 + Y.^2) < 15;

fat = im;
fat(180:190,100:200) = 1;
%}
load('pd_phantom.mat')
pd = rot90(pd,-1);
heart = rot90(heart,-1);

phantom1 = pd;
phantom2 = ~heart.*pd2 + heart.*pd;

para.angle_mod = 1;
para.Recon.sy = 20;
para.Recon.sz = 1;
para.dataAngle = 180;

theta = get_angle_mod(para);
[kx,ky] = get_k_coor(288,theta,0,144);
N = NUFFT.init_new(kx,ky,1.5,[6 6]);

kSpace = NUFFT.NUFFT_new(phantom1,N);
im_recon = NUFFT.NUFFT_adj_new(kSpace,N);

kSpace2 = NUFFT.NUFFT_new(phantom2,N);
im_recon2 = NUFFT.NUFFT_adj_new(kSpace2,N);

sig1 = sum(sum(abs(heart .* im_recon)))
sig1 = sum(sum(abs(heart .* im_recon2)))

figure,subplot(2,3,1)
imagesc(phantom1)
colormap gray
brighten(0.4)
axis image
axis off
title 'test image'

subplot(2,3,2)
imagesc(abs(im_recon))
colormap gray
brighten(0.4)
axis image
axis off
title '#1 recon'

subplot(2,3,3)
imagesc(abs(im_recon2))
colormap gray
brighten(0.4)
axis image
axis off
title '#2 recon with increased fat signal'

subplot(2,3,4)
imagesc(abs(heart.*pd))
colormap gray
brighten(0.4)
axis image
axis off
title 'blood pool mask'

subplot(2,3,5)
imagesc(abs(fat.*pd))
colormap gray
brighten(0.4)
axis image
title 'fat mask'
axis off