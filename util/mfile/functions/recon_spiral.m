function recon_spiral(data_dir)

NFOV = 2.5;
weight_tTV = 0.08 / 4;
weight_sTV = 0.00;
nos_one = 1; % number of spirals per time frame
%%
imsize = [84, 84] * NFOV;

load(fullfile(data_dir.folder, data_dir.name));

scale_factor = 1e3*prod(imsize)/max(abs(kdata(:)));

kSpace = single(permute(kdata,[1,3,2]))*scale_factor;

clear kdata
clear spokeindex

sx = size(kSpace,1);
nos = size(kSpace,2);
nc = size(kSpace,3);

kx = real(kloc)*imsize(1);
ky = imag(kloc)*imsize(2);

clear kloc

Nframe = size(kSpace,2)/nos_one;
Nframe = floor(Nframe);
% Nframe = 100; % hard coded to reconstruct only 30 time frames for testing.
frame_begin = 1;
% drop the data at the end
kSpace = kSpace(:,(frame_begin-1)*nos_one + 1 : (frame_begin + Nframe - 1)*nos_one,:);
kx = kx(:,(frame_begin-1)*nos_one + 1 : (frame_begin + Nframe - 1)*nos_one);
ky = ky(:,(frame_begin-1)*nos_one + 1 : (frame_begin + Nframe - 1)*nos_one);
nos = Nframe * nos_one;
w = repmat(w,[1,ceil(nos/size(w,2))]);
w(:,Nframe*nos_one+1:end) = [];
kx(:,Nframe*nos_one+1:end) = [];
ky(:,Nframe*nos_one+1:end) = [];
nos = Nframe*nos_one;

%% NUFFT
kSpace = reshape(kSpace,[sx,nos_one,Nframe,nc]);
kx = reshape(kx,[sx,nos_one,Nframe]);
ky = reshape(ky,[sx,nos_one,Nframe]);
N = NUFFT.init(kx,ky,1,[4,4],imsize(1),imsize(1));
N.W = single(w(:,1));

Image = NUFFT.NUFFT_adj(kSpace,N);

sens = get_sens_map(Image,'2D');

Data.kSpace = kSpace;
Data.N = N;
Data.sens_map = sens;
Data.first_est = Image.*conj(sens);
Data.first_est = sum(Data.first_est,4);
scale = max(abs(Image(:)));
clear Image sens kSpace kx ky w

%% Recnstruction parameters
para.setting.ifplot = 1;
para.setting.ifGPU = 1;
para.Recon.epsilon = eps('single');
para.Recon.step_size = 2;
para.Recon.ifContinue = 0;
para.Recon.noi = 150; % number of iterations
para.Recon.type = '2D Spiral server'; %stack of spiral
para.Recon.no_comp = nc;
para.Recon.break = 1;

para.Recon.weight_tTV = scale * weight_tTV; % temporal regularization weight
para.Recon.weight_sTV = scale * weight_sTV; % spatial regularization weight

% clearvars -except Data para iraw weight_tTV

%% conjugate gradient solver
[Image_recon_ar, para] = STCR_conjugate_gradient(Data,para);

Image_recon_ar = abs(rot90(crop_half_FOV(Image_recon_ar, [84, 84])));

save_dir = ['./ReconData/', data_dir.folder(end-10:end), '/'];
if isempty(dir(save_dir))
    mkdir(save_dir)
end

save_dir = [save_dir, data_dir.name];

save(save_dir, 'Image_recon_ar', 'para')
end