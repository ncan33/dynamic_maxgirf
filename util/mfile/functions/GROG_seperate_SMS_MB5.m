function [G,kSpace_cart,kSpace_radial_out] = GROG_seperate_SMS_MB5(kSpace_radial, kx, ky, phase_mod, para)
%  This is a GROG interpolation function for 2D MRI k-space data

%  Reference:
%  Nicole Seiberlich, et al (2008) Magnetic Resonance in Medicine
%  Self-Calibrating GRAPPA Operator Gridding for Radial and Spiral
%  Trajectories. 59:930-935.

%  Inputs:
%         kSpace_radial: up to 5 dimensions
%         kx, ky       : k-space positions

%  Copyright. Ye Tian, U of U DiBella Group
%  phye1988@gmail.com


% This version does not support multi-slices inputs.
nSMS = para.Recon.nSMS;
[sx,nor,nof,nc,~,NSlice] = size(kSpace_radial);
core_size = para.core_size;
over_sampling = para.over_sampling;

if nSMS == 5
    phase_mod = phase_mod(:,:,:,:,2);
end

phase_mod = mod(angle(phase_mod)/(2*pi/5),nSMS);

SMS(1,:,:,1,1,1) = phase_mod == 0;
SMS(1,:,:,1,1,2) = phase_mod == 1;
SMS(1,:,:,1,1,3) = phase_mod == 2;
SMS(1,:,:,1,1,4) = phase_mod == 3;
SMS(1,:,:,1,1,5) = phase_mod == 4;

G = cell(1,nSMS);
for j=1:nSMS
    for i=1:NSlice
        kSpace_radial_temp = kSpace_radial(:,:,:,:,i);
        kSpace_radial_temp = kSpace_radial_temp(repmat(SMS(:,:,:,:,:,j),[sx 1 1 nc]));
        kx_temp = kx(repmat(SMS(:,:,:,:,:,j),[sx 1 1]));
        ky_temp = ky(repmat(SMS(:,:,:,:,:,j),[sx 1 1]));
        kSpace_radial_temp = reshape(kSpace_radial_temp,[sx nor/nSMS nof nc]);
        kx_temp = reshape(kx_temp,[sx nor/nSMS nof]);
        ky_temp = reshape(ky_temp,[sx nor/nSMS nof]);
        G{j} = GROG.GNUFFT_init(kSpace_radial_temp,kx_temp,ky_temp,over_sampling,core_size);
        kSpace_cart(:,:,:,:,:,:,j) = GROG.GNUFFT_rad2cart(kSpace_radial_temp,G{j});
        kSpace_radial_out(:,:,:,:,:,:,j) = kSpace_radial_temp;

    end
end

kSpace_radial_out = reshape(kSpace_radial_out,[sx,nor/nSMS,nof,nc,1,NSlice,nSMS]);
sx_over = size(kSpace_cart,1);
kSpace_cart = reshape(kSpace_cart,[sx_over,sx_over,nof,nc,1,NSlice,nSMS]);
