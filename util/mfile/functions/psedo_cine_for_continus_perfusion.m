para.Recon.sy=340;
para.Recon.sz=70;
para.dataAngle=180;
para.angle_mod = 1;
para.phase_mod = 1;

theta = get_angle_mod(para);
theta = reshape(theta,[1,340,70]);
phase_mod = get_phase_mod(para);
phase_mod = reshape(phase_mod,1,340,70,3);
[kx,ky] = get_k_coor(288,theta,0,144);


nof = 100;
window_length = 36;
sliding_length = 3;

SW = 1:window_length;
SW = repmat(SW,nof,1);
SW = bsxfun(@plus,SW,(0:sliding_length:sliding_length*(nof-1)).');

im = zeros(288,288,nof,8);

for i=1:nof
tic
frame1 = kSpace(:,SW(i,:),:,:);
frame1 = reshape(frame1,288,70*window_length,1,8);
frame1_kx = kx(:,SW(i,:),:); frame1_kx = reshape(frame1_kx,288,70*window_length,1);
frame1_ky = ky(:,SW(i,:),:); frame1_ky = reshape(frame1_ky,288,70*window_length,1);
N = NUFFT.init(frame1_kx,frame1_ky,1.5,[6 6]);
im(:,:,i,:) = NUFFT.NUFFT_adj(frame1,N);
toc
end
im_for_sens = squeeze(sum(im,3));
sens_map = ismrm_estimate_csm_walsh_optimized_yt(im_for_sens);
sens_map = permute(sens_map,[1 2 4 3]);
im_sens = bsxfun(@times,im,conj(sens_map));
im_sens = sum(im_sens,4);

%im_sos = conj(im).*im;
%im_sos = sum(im_sos,4);
%im_sos = sqrt(im_sos);