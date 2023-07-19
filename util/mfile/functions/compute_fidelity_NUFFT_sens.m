function fidelity_update = compute_fidelity_NUFFT_sens(im,kSpace_radial,N_f)

k = fft2(im.*N_f{1}.Apodizer);

for i=1:size(kSpace_radial,4)
    k_radial = NUFFT.cart2rad_new(k,N_f{i})/10^5;
    d = kSpace_radial(:,:,:,i) - k_radial;
    d(1:72,:,:) = 0;
    k_cart(:,:,:,i) = NUFFT.rad2cart(d.*N_f{i}.W,N_f{i});
end
k_cart = sum(k_cart,4)/4.5619e+04;
fidelity_update = ifft2(k_cart);
end