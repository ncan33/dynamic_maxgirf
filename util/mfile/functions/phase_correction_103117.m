function kSpace = phase_correction_103117(kSpace)

nor = size(kSpace,2);
mod0 = 1:3:nor;
center = kSpace(145,mod0,:,1,:);
center_0 = kSpace(145,:,:,1,:); % save for plot

center_phase_frame_1 = angle(center(1,:,1,1,1));
center_mean_frame_1 = mean(center_phase_frame_1,2);
center_phase_frame_1 = center_phase_frame_1 - center_mean_frame_1;


phase_all = -pi:0.01:pi;
phase_all = exp(1i*phase_all).';

phase_shift = phase_all.*(center);
phase_shift_diff = angle(phase_shift) - center_phase_frame_1;
phase_shift_diff_sos = sum(phase_shift_diff.^2,2);
[~,idx] = min(phase_shift_diff_sos);
phase_correction = phase_all(idx);
if size(phase_correction,1)~=1
    phase_correction = permute(phase_correction,[2 3 1]);
end
kSpace = phase_correction.*kSpace;


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