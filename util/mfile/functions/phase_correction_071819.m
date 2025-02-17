function kSpace = phase_correction_071819(kSpace)
[nx,nor,nof,nc,nset] = size(kSpace);
center_idx = round(nx/2)+1;

coil_idx = squeeze(sum(sum(sum(abs(kSpace(center_idx,:,:,:,:)).^2,2),3),5));
[~,coil_idx] = max(coil_idx);

%kSpace_temp = reshape(kSpace,nx,nor*nof,nc,nset);
% mod0 = phase_mod == 0;
%mod0 = false(1,nor*nof);
%mod0(1:3:nor*nof) = true;
% mod0 = reshape(mod0,1,nor,nof,1,size(mod0,2));
% N = floor(nor/3)*3;
% mod0(:,N+1:end,:,:,:) = false;

center = kSpace(center_idx,:,:,coil_idx,:);
%center = reshape(center,1,N/3,nof,1,nset);
center_0 = center; % save for plot
% keyboard
%center_phase_frame_1 = angle(center(1,:,1,1,1));
%center_phase_frame_1 = mod(center_phase_frame_1,2*pi);
%center_mean_frame_1 = median(center_phase_frame_1);
%center_phase_frame_1 = center_phase_frame_1 - center_mean_frame_1;

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


center = kSpace(center_idx,:,:,coil_idx,:);
figure,plot(angle(center_0(:)))
hold on
plot(angle(center(:)))
title 'Phase Correction Result'
xlabel 'phase mod 0 measurements'
ylabel 'phase'
legend 'before' 'after'
drawnow