function [G,kSpace_cart,kSpace_radial_out] = GROG_seperate_SMS_GNUFFT_MB2( kSpace_radial, kx, ky, phase_mod, para)
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


phase_mod = phase_mod(:,:,:,:,2);


idx = phase_mod==0;
%phase_mod = round(imag(phase_mod));
% temp = [0,1,-1];
% N = sum(idx,2)/3;
% for i = 1:size(idx,3)
%     phase_mod(1,idx(1,:,i),i) = repmat(temp,[1,N(i)]);
% end
SMS(1,:,:,1,1,1) = phase_mod==1;
SMS(1,:,:,1,1,2) = phase_mod==-1;
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
        
        %[kx_temp,ky_temp] = GROG.Trajectory_correction_new(kSpace_radial_temp,kx_temp,ky_temp,400);
        G{j} = GROG.GNUFFT_init(kSpace_radial_temp,kx_temp,ky_temp,over_sampling,core_size);
        kSpace_cart(:,:,:,:,:,:,j) = GROG.GNUFFT_rad2cart(kSpace_radial_temp,G{j});
        kSpace_radial_out(:,:,:,:,:,:,j) = kSpace_radial_temp;
    
    
    %if np==4
  
        %[kx_temp1,ky_temp1] = GROG.Trajectory_correction_DC_sum(kSpace_radial_temp,kx_temp,ky_temp,50);

        %

        
        %G{j} = GROG.GNUFFT_init(kSpace_radial_temp,kx_temp,ky_temp,over_sampling,core_size);

        

        %kSpace_cart(:,:,:,:,:,:,j) = GROG.GNUFFT_rad2cart(kSpace_radial_temp.*G{j}.W,G{j});
        %mask(:,:,:,:,:,:,j) = GROG.GNUFFT_rad2cart(ones(size(kSpace_radial_temp)).*G{j}.W,G{j});

        %N{j} = NN.NN_init(kx_temp,ky_temp,para.over_sampling);
        %mask(:,:,:,:,:,:,j) = NN.NN_rad2cart(ones(size(kSpace_radial_temp)).*G{j}.W,N{j});
        
    %elseif np==1
    %    G{j} = GROG.GROG_init(kSpace_radial_temp,kx_temp,ky_temp);
    %    kSpace_cart(:,:,:,:,:,:,j) = GROG.GROG_rad2cart(kSpace_radial_temp,G{j});
    %end

    %{
    This sucks! The GROG operator must be super non-perfect! Learn the GROG
    trajectory correction!
    
    noi = 100;
    step = 0.1;
    k_cart   = GROG.GROG_rad2cart(kSpace_radial_temp,G{1});

for i=1:noi
    
    k_radial = GROG.GROG_cart2rad(k_cart,G{1});
    k_radial_updte = kSpace_radial_temp - k_radial;
    k_cart_update = GROG.GROG_rad2cart(k_radial_updte,G{1});
    k_cart = k_cart + step*k_cart_update;
    
    norm = abs(k_cart_update);
    norm = sqrt(sum(norm(:).^2));
    norm_all(i) = norm;
    figure(1)
    plot(norm_all)
    %figure(2)
    %imagesc(abs(k_cart(:,:,1,1)))
    drawnow
end
    kSpace_cart(:,:,:,:,:,:,j) = k_cart;
    %}
    end
end
%kSpace_radial_out = mask;
kSpace_radial_out = reshape(kSpace_radial_out,[sx,nor/nSMS,nof,nc,1,NSlice,nSMS]);
sx_over = size(kSpace_cart,1);
kSpace_cart = reshape(kSpace_cart,[sx_over,sx_over,nof,nc,1,NSlice,nSMS]);
kSpace_cart([1,end],:,:,:,:,:,:) = [];
kSpace_cart(:,[1,end],:,:,:,:,:) = [];
