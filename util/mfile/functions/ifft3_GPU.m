function imspace = ifft3_GPU(kspace)
kspace = gpuArray(kspace);
imspace = ifft(ifft(ifft(kspace,[],1),[],2),[],3);
%imspace = fftshift(imspace,3);
imspace = gather(imspace);
end