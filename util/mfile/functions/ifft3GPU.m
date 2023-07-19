function imspace = ifft3GPU(kspace)

imspace(:,:,:,:,1:4) = ifft2(kspace(:,:,:,:,1:4));
imspace(:,:,:,:,1:4) = ifft(imspace(:,:,:,:,1:4),[],3);
imspace(:,:,:,:,5:8) = ifft2(kspace(:,:,:,:,5:8));
imspace(:,:,:,:,5:8) = ifft(imspace(:,:,:,:,5:8),[],3);
%imspace = ifft(ifft(ifft(kspace,[],1),[],2),[],3);
%imspace = fftshift(imspace,3);

end