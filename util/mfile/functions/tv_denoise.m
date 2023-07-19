function im_denoised = tv_denoise(im,weight,noi)
im_denoised = im;
%figure(1)
for i=1:noi
    fprintf(sprintf('iter = %d \n',i))
    im_denoised = im_denoised + 0.1 * compute_sTV_yt(im_denoised,weight,1e-6) + 0.1 * (im-im_denoised);
    %imagesc(abs(im_denoised(:,:,1)))
    %colormap gray
    %axis image
    %axis off
    %brighten(0.4)
    %drawnow
end