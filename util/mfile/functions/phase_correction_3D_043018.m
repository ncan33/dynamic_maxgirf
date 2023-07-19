function kSpace = phase_correction_3D_043018(kSpace,slice_idx,measurements)
[nx,~,~,~] = size(kSpace);
center_idx = round(nx/2)+1;

coil_idx = squeeze(sum(abs(kSpace(center_idx,:,:,:,:)).^2,2));
[~,coil_idx] = max(coil_idx);

%PD = ray_mod == 0;

%kSpace_temp = reshape(kSpace,nx,nor*nof,nc,nset);
if min(slice_idx(:)) <= 0
    slice_idx = slice_idx - min(slice_idx) + 1;
end
kz0 = slice_idx == 5;
%mod0 = false(1,nor*nof);
%mod0(1:3:nor*nof) = true;
%kz0 = reshape(kz0,1,nor,1,1);
%N = floor(nor/3)*3;
%kz0(:,N+1:end,:,:,:) = false;

center = kSpace(center_idx,:,:,coil_idx);
phase_all = -pi:0.01:pi;
phase_all = exp(1i*phase_all).';

for i=1:max(measurements)
    idx_measurement = measurements == i;
    center_temp = center(idx_measurement & kz0);
    phase_shift = phase_all.*(center_temp);
    phase_shift_diff = angle(phase_shift);
    phase_shift_diff_sos = sum(phase_shift_diff.^2,2);
    [~,phase_shift_idx] = min(phase_shift_diff_sos);
    phase_correction = phase_all(phase_shift_idx);
    kSpace(:,idx_measurement,:,:) = phase_correction.*kSpace(:,idx_measurement,:,:);
end
%center = reshape(center,1,N/3,nof,1,nset);
%center_0 = kSpace(center_idx,:,:,coilidx,:); % save for plot

%center_phase_frame_1 = angle(center(1,:,1,1,1));
%center_phase_frame_1 = mod(center_phase_frame_1,2*pi);
%center_mean_frame_1 = median(center_phase_frame_1);
%center_phase_frame_1 = center_phase_frame_1 - center_mean_frame_1;



% - center_phase_frame_1;


%phase_correction;
%if size(phase_correction,1)~=1
%    phase_correction = permute(phase_correction,[2 3 1]);
%end
%kSpace = phase_correction.*kSpace;


%center = kSpace(145,:,:,1,:);
%figure,plot(angle(center(:)))
%hold on
%plot(angle(center_0(:)))
%drawnow










%phase_shift = mean(angle(phase_shift),2);
%[~,idx] = min(abs(phase_shift));
%phase_correction = phase_all(idx);

%kSpace = phase_correction.*kSpace;

%center = kSpace(145,:,:,1,:);
%center = center(:,:,:,:,1);
%figure,plot(angle(center(:)))
