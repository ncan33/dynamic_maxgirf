function kSpace = phase_correction_022718(kSpace)
[nx,nor,nof,nc,nset] = size(kSpace);
center_idx = round(nx/2)+1;
coilidx = squeeze(sum(sum(sum(abs(kSpace(center_idx,:,:,:,:)).^2,2),3),5));
[~,coilidx] = max(coilidx);


kSpace_temp = reshape(kSpace,nx,nor*nof,nc,nset);
mod0 = false(1,nor*nof);
mod0(1:3:nor*nof) = true;
mod0 = reshape(mod0,nor,nof);
N = floor(nor/3)*3;
mod0(N+1:end,:) = false;

center = kSpace_temp(center_idx,mod0,coilidx,:);
center = reshape(center,1,N/3,nof,1,nset);
%center_0 = kSpace(center_idx,:,:,coilidx,:); % save for plot

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
