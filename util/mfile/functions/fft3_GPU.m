function kspace = fft3_GPU(imspace)
imspace = gpuArray(imspace);
%imspace = fftshift(imspace,3);
kspace = fft(fft(fft(imspace,[],1),[],2),[],3);
kspace = gather(kspace);
end