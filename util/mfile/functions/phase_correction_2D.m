function kSpace = phase_correction_2D(kSpace)
[nx,nor,nof,nc,nset] = size(kSpace);
center_idx = round(nx/2)+1;

coil_idx = squeeze(sum(sum(abs(kSpace(center_idx,:,:,:,:)).^2,2),3));
[~,coil_idx] = max(coil_idx);

center = kSpace(center_idx,:,:,:,coil_idx);
center_0 = center; % save for plot

phase_all = -pi:0.01:pi;
phase_all = exp(1i*phase_all).';

phase_shift = phase_all.*(center);
phase_shift_diff = angle(phase_shift);% - center_phase_frame_1;
phase_shift_diff_sos = sum(phase_shift_diff.^2,2);
[~,idx] = min(phase_shift_diff_sos);
phase_correction = phase_all(idx);
if size(phase_correction,1)~=1
    phase_correction = permute(phase_correction,[2 3 1]);
end
kSpace = phase_correction.*kSpace;


center = kSpace(center_idx,:,:,:,coil_idx);
figure,plot(angle(center_0(:)))
hold on
plot(angle(center(:)))
title 'Phase Correction Result'
xlabel 'phase mod 0 measurements'
ylabel 'phase'
legend 'before' 'after'
drawnow










%phase_shift = mean(angle(phase_shift),2);
%[~,idx] = min(abs(phase_shift));
%phase_correction = phase_all(idx);

%kSpace = phase_correction.*kSpace;

%center = kSpace(145,:,:,1,:);
%center = center(:,:,:,:,1);
%figure,plot(angle(center(:)))
