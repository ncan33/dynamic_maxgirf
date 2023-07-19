function Image = get_cardiac_phase_continues_SMS_all_slices(file)
load(fullfile(file.folder,file.name))

set = update_SET;
theta = kSpace_info.angle_mod;
phase = kSpace_info.phase_mod;
frames = kSpace_info.frames;
PD_frames = kSpace_info.ProtonDensityScans;

nof = max(frames(:));
sx = size(kSpace,1);
nc = size(kSpace,5);
nor = diff(find([1,diff(frames),length(theta)]));
kSpace_radial = zeros(sx,1,nof-PD_frames,nc,'single');
theta_all = zeros(1,1,nof-PD_frames,'single');
phase_all = repmat([0,1,2],[1,floor(max(nor)/3),nof-PD_frames]);

for i = PD_frames+1:nof
    rays = (1:nor(i)) + sum(nor(1:i-1));
    kSpace_temp = squeeze(kSpace(:,rays,:,:,:));
    theta_temp = theta(rays);
    phase_temp = phase(rays);
    rays_idx = 1:nor(i);
    nor_temp = floor(rays_idx(end)/3);
    rays_idx(nor_temp*3+1:end) = [];
    kSpace_radial(:,rays_idx,i-PD_frames,:) = kSpace_temp(:,rays_idx,:);
    theta_all(:,rays_idx,i-PD_frames) = theta_temp(:,rays_idx,:);
    phase_all(:,rays_idx,i-PD_frames) = phase_temp(:,rays_idx,:);
end

window_length = 9;

nof_sl = min(floor(median(nor)/9),100);
Image = zeros(sx,sx,nof_sl,nc,3,'single');

for i=1:nof_sl
    rays_idx = (1:window_length) + (i-1)*window_length;
    kSpace_temp = kSpace_radial(:,rays_idx,:,:);
    kSpace_temp = reshape(kSpace_temp,[sx,window_length*(nof-PD_frames),1,1,nc]);
    theta_temp = theta_all(1,rays_idx,:);
    phase_temp = phase_all(1,rays_idx,:);
    phase_temp = phase_temp(:).';
    clear phase_mod_temp
    phase_mod_temp(1,:,1,1,1,2) = exp(1i*2*pi/3*phase_temp);
    phase_mod_temp(1,:,1,1,1,3) = exp(1i*4*pi/3*phase_temp);
    phase_mod_temp(1,:,1,1,1,1) = 1;
    [kx,ky] = get_k_coor(sx,theta_temp(:).',0,round(sx/2+1));
    N = NUFFT.init_new(kx,ky,1.5,[6,6]);
    Image(:,:,i,:,:) = NUFFT.NUFFT_adj_new(kSpace_temp.*conj(phase_mod_temp),N);
end

Image = squeeze(crop_half_FOV(abs(sos(Image,4))));
MID = strfind(file.name,'MID');
MID = file.name(MID+5:MID+7);
folder = ['ReconData/mat_files/MID',MID,'/'];
mkdir(folder)
save([folder,'SET',num2str(set),'_',file.name(6:end-4),'_cine_sliding_9.mat'],'Image')


% window_length = 40;
% Image_sw = zeros(sx,sx,nof_sl-window_length+1,nc,'single');
% 
% for i=1:nof_sl-window_length+1
%     rays_idx = (1:window_length) + i-1;
%     Image_sw(:,:,i,:) = sum(Image(:,:,rays_idx,:),3);
% end
% 
% Image_sw = crop_half_FOV(sos(Image_sw));
% 
% keyboard