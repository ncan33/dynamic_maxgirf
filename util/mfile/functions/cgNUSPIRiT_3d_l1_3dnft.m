function [x_old, norm, s_norm, f_norm] = cgNUSPIRiT_3d_l1_3dnft(kspace, x0, N, kernel, niter, lambda, tv_weight)

% Implementation of image-domain SPIRiT reconstruction from arbitrary
% k-space. The function is based on Jeff Fessler's nufft code and LSQR 
% 
% Inputs: 
%       kData   - k-space data matrix it is 3D corresponding to [readout,interleaves,coils] 
%       x0      - Initial estimate of the coil images
%       NUFFTOP - nufft operator (see @NUFFT class)
%       GOP     - SPIRiT Operator (See @SPIRiT)
%       nIter   - number of LSQR iterations
%       lambda  - ratio between data consistency and SPIRiT consistency (1 is recommended)
%
% Outputs:
%       res - reconstructed coil images
%       FLAG,RELRES,ITER,RESVEC,LSVEC - See LSQR documentation
%
% See demo_nuSPIRiT for a demo on how to use this function.
%
% (c) Michael Lustig 2006, modified 2010
%% 
N.W     = sqrt(N.W);
kspace  = kspace .* N.W;
kspace  = permute(kspace, [1, 2, 4, 3]);
sx      = size(x0, 1);
sy      = size(x0, 2);
sz      = size(x0, 3);
ncoil   = size(x0, 4);
x0      = x0 * sqrt(sx * sy * sz);

% undersample_mask = logical(abs(kspace(1, :, :, :, 1)));

%%
imSize      = size(x0);
% dataSize    = size(kspace);
x_old       = x0;
step        = 0.35;
tv_weight   = max(abs(x0(:))) * tv_weight;

x_old   = gpuArray(x_old);
N.S     = gpuArray(N.S);
% y1      = gpuArray(zeros(dataSize, 'single'));
% kernel  = gpuArray(single(kernel));

s_norm = zeros([1, niter]);
f_norm = zeros([1, niter]);
norm   = zeros([1, niter]);

for i = 1:niter
    tic
    %% spirit update
    
    s_update = squeeze(sum(kernel .* x_old, 4)) - x_old;
    
    s_norm(i)= sum(abs(s_update(:)).^2) * lambda;
    
    s_update = sum(conj(kernel) .* permute(s_update, [1, 2, 3, 5, 4]), 5) - s_update;
    
%     s_update = 0;
%     s_norm(i) = 0;
    %% fidelity update
    f_update = NUFFT3.NUFFT(x_old, N);
    f_update = f_update / sqrt(prod(imSize(1:3))) .* (N.W);
    
    f_update = kspace - f_update;
    
    f_norm(i)= sum(abs(f_update(:)).^2);
    
    f_update = NUFFT3.NUFFT_adj(f_update, N);
    f_update = squeeze(f_update) * sqrt(prod(imSize(1:3)));
    
    %% tv update
    tv_update = zeros(size(x_old), 'like', x_old);
    for icoil = 1:ncoil
        tv_update(:, :, :, icoil) = compute_sTV_3D_yt(x_old(:, :, :, icoil), tv_weight, eps('single'));
    end
    %% cg update
%     x_new = x_old - s_update * 0.1 + f_update * 0.1;
%     x_old = x_new;
    norm(i) = f_norm(i) + s_norm(i);
    
    update_term = (f_update * 0.2 - s_update * lambda + tv_update);
    
    if i > 1
        beta = update_term(:)'*update_term(:)/(update_term_old(:)'*update_term_old(:)+eps('single'));
        update_term = update_term + beta*update_term_old;
        if norm(i) > norm(i-1)
            step = step / 2;
            step
        end
    end
    
    x_old = x_old + update_term * step;
    update_term_old = update_term;
    
    
    figure(1)
    imagesc(abs(sos(x_old(:,:,round(sx/2), :))))
    axis image
    title(i)
    drawnow
    
    figure(2)
    clf
    loglog(norm)
    hold on
    plot(f_norm)
    plot(s_norm)

    toc
end
x_old = gather(x_old);