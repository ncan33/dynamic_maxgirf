function im_denoised = ttv_denoise(im,weight,noi)
im_denoised = im;
figure(1)
for i=1:noi
    %im_denoised = im_denoised + 0.1 * compute_sTV_yt(im_denoised,weight,1e-6);
    im_denoised = im_denoised + 0.1 * compute_tTV_yt(im_denoised,weight,1e-6);
    imagesc(im_denoised(:,:,50))
    colormap gray
    axis image
    axis off
    brighten(0.4)
    drawnow
end