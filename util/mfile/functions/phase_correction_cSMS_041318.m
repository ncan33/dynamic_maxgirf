function kSpace = phase_correction_cSMS_041318(kSpace,kSpace_info)

frames = kSpace_info.frames;
phase_mod = kSpace_info.phase_mod;
idx_mod0 = phase_mod == 0;
nof = max(frames);
center = round(size(kSpace,1)/2+1);

phase_all = -pi:0.01:pi;
phase_all = exp(1i*phase_all).';

correction = zeros(1,length(frames),'single');


[~,coil_idx] = max(abs(squeeze(sum(kSpace(center,:,:)))));

for i=1:nof
    idx_frames = frames == i;
    idx_temp = idx_frames & idx_mod0;
    phase_temp = kSpace(center,idx_temp,coil_idx);
    phase_temp = phase_temp(1:min(10,length(phase_temp)));
    phase_shift = phase_temp.*phase_all;
    phase_shift = angle(phase_shift);
    phase_shift = sum(phase_shift.^2,2);
    [~,idx_phase] = min(phase_shift);
    correction(idx_frames) = phase_all(idx_phase);
end

kSpace = kSpace.*correction;