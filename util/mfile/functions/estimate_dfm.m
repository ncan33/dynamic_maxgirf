function Data =  estimate_dfm(Data, para)
%   [dfm, mask, coil_comp_img] =  estimate_dfm(data, csm, para)
%
%   Estimates dynamic field map from signle echo-time spiral data
%
%   INPUT:
%     - data.w      [samples, full_spirals]   : density compensation func
%     - data.kdata  [samples, spirals, coil]  : spiral k-space 
%     - data.kloc   [samples, spirals]        : k-space coordinate
%     - csm         [x,y,coil]                : Relative coil sensitivity maps
%     - para
%
%   OUTPUT:
%     - dfm         [x,y,frames]    : dynamic field map estimate
%     - mask        [x,y,frames]    : binary mask
%     - c_comp_img  [x,y,frames]    : coil composited dynamic image
%
% (c) Yongwan Lim (yongwanl@usc.edu) University of Southern California, 2018.
%
% Modified Ye Tian (phye1988@gmail.com) 03/08/2021

crop_thresh = para.crop_thresh; % 1/para.te/2;

%% sliding-window for dfm estimation
[sx, narm, nof, nc] = size(Data.kSpace);
Data.kx = reshape(Data.kx, [sx, 1, narm * nof]);
Data.ky = reshape(Data.ky, [sx, 1, narm * nof]);
kSpace = reshape(Data.kSpace, [sx, 1, narm * nof, nc]);
N = NUFFT.init(Data.kx, Data.ky, 1, [4, 4], para.Recon.image_size(1), para.Recon.image_size(1));
N.W = Data.N.W;
Image_one_arm = NUFFT.NUFFT_adj(kSpace, N);

Image_sw = zeros([para.Recon.image_size, nof, nc], 'single');

arm_off_idx = round((para.NarmsFull - para.Recon.narm) / 2);

for iframe = 1:nof
    idx = (1 : para.NarmsFull) + (iframe - 1) * para.Recon.narm - arm_off_idx;
    if idx(1) < 1
        idx = 1 : para.NarmsFull;
    end
    if idx(end) > narm * nof
        idx = narm * nof + 1 - para.NarmsFull: narm * nof;
    end
    Image_sw(:, :, iframe, :) = sum(Image_one_arm(:, :, idx, :), 3);
end

%% coil combine
% the coil estimation is different than the coil estimation used in
% reconstruction. This is based on the SOS coil combination, and offers a
% phase-removed image for off-resonance map estimaiton.

kernel = gen_2d_kernel(para.sizehanningfilter, 2); % spatial filter window (Hanning window)
normalization = conv2(ones(para.Recon.image_size), kernel,'same'); % normalization mask

talr_img = imfilter(sum(Image_one_arm, 3), kernel) ./ normalization;
% sens_map = ismrm_estimate_csm_mckenzie(squeeze(talr_img), 0);
sens_map = ismrm_estimate_csm_mckenzie(squeeze(talr_img));
sens_map = permute(sens_map, [1, 2, 4, 3]);
Image_sw = sum(Image_sw .* conj(sens_map), 4);

%% Smooth the coil-composited image

% spatial smoothing
kernel_s = gen_2d_kernel(3, 2);
normalization_mask = 1 ./ convn(ones(para.Recon.image_size),kernel_s,'same'); % normalization mask
c_comp_img_s = imfilter(Image_sw, kernel_s, 'same') .* normalization_mask;

% temporal smoothing
c_comp_img_s_t = circshift(c_comp_img_s, [0, 0, 1]) * 0.25 + circshift(c_comp_img_s, [0, 0, -1]) * 0.25 + c_comp_img_s * 0.5;
% c_comp_img_s_t = c_comp_img_s; % YT: 2021/05/28

%% Estimate dynamic field map
dfm = angle(c_comp_img_s_t) / (-2 * pi * para.te);

dfm(isnan(dfm)) = 0;
dfm(dfm >=  crop_thresh) = crop_thresh;
dfm(dfm <= -crop_thresh) = -crop_thresh;

%% Mask for field map 
mask = zeros(size(dfm));
thresh = 0.01;
for tt = 1:nof
    temp = abs(c_comp_img_s_t(:,:,tt)).^2;
    temp = temp/max(temp(:));    
    temp(temp<thresh) = 0;
    temp(temp>thresh) = 1;
    se = strel('disk',1);
    temp = imclose(temp,se);
    mask(:, :, tt) = temp;
end

c_comp_img = Image_sw./max(abs(Image_sw(:)));

Data.off_res.map = dfm;
% Data.off_res.map = dfm .* mask; % YT: 2021/05/28
Data.sens_map = sens_map;
