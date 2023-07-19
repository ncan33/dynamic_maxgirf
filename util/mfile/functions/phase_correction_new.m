function kSpace = phase_correction_new(kSpace)

phase_corr = (-pi:0.01:pi).';
phase_corr = exp(1i*phase_corr);
center_line = kSpace(145,:,:,:,:);
[~,coil_selected] = min(std(angle(squeeze(center_line(:,:,:,1,1)))));
center_line = center_line(:,:,coil_selected,:,:);

center_line_phase_corr = center_line.*phase_corr;
center_line_phase_avrg = mean(angle(center_line_phase_corr),2);
%center_line_phase_avrg(center_line_phase_avrg<0) = inf;
[~,phase_position] = min(abs(center_line_phase_avrg));
%phase_position(phase>1);
phase_corr = phase_corr(phase_position);
kSpace = kSpace.*phase_corr;
center_line_2 = kSpace(145,:,:,:,:);
center_line_phase_max = max(angle(center_line_2));
pi_corr = center_line_phase_max>2.5;
kSpace = kSpace.*exp(-1i*pi_corr*pi);



center_line_before = squeeze(center_line);
center_line_after  = squeeze(kSpace(145,:,coil_selected,:,:));
%center_line_before = permute(center_line_before,[1 3 2 4]);
%center_line_after  = permute(center_line_after,[1 3 2 4]);
siz = size(center_line_after);
center_line_before = reshape(center_line_before,[siz(1)*siz(2),siz(3)]);
center_line_after  = reshape(center_line_after,[siz(1)*siz(2),siz(3)]);

figure
plot(angle(center_line_before(:,1,1)));
hold on
plot(angle(center_line_after(:,1,1)));
