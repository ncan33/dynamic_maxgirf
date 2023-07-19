function vbm4d_denoised = VBM4D_ye(image)

%realImg(:,:,:,1) = real(image);
%realImg(:,:,:,2) = imag(image);
abs_img = abs(image);
phase = angle(image);
for i=1:size(image,4)
    denoisedImg(:,:,:,i) = VBM4D.vbm4d(abs_img(:,:,:,i));
end

vbm4d_denoised = denoisedImg + exp(1i*phase);
%vbm4d_denoised = permute(vbm4d_denoised,[1 2 3 5 4]);


return;