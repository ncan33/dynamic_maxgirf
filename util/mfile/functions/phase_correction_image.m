function kSpace = phase_correction_image(image, kSpace)
[sx, sy, nof, ncoil, nslice] = size(image);
phase = -pi:0.01:pi;
phase = exp(1i*phase);
phase = reshape(phase,1,1,629);

% phase = gpuArray(phase);
%image = gpuArray(image);
phase_mean = angle(sum(image(:,:,:,1,:), 3));
for i=1:nof
    for j=1
        for k=1:nslice
%             phase_standard = angle(image(:,:,i,j,k));
            image_all_phase = image(:,:,i,j,k).*phase;
            phase_all = angle(image_all_phase);
            phase_difference = phase_all - phase_mean;
            phase_difference_sum = sum(sum(abs(phase_difference)));
            [~,loc] = min(phase_difference_sum,[],3);
            phase_correction(1,1,i,1,k) = phase(loc);
%             image(:,:,i,j,k) = gather(image_all_phase(:,:,loc));
        end
    end
end
kSpace = kSpace .* phase_correction;
%image = gather(image);
