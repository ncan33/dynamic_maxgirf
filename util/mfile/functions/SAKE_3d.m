function im = SAKE_3d(kspace, N, kernel_size, wnRank, nIter, show)
%
% res = SAKE(DATA, kSize, wnRank,nIter)
%
% function performs an autocalibration reconstrution without calibration
% lines based on the SAKE method (P. Shin et. al, "Calibrationless Parallel 
% Imaging Reconstruction Based on Structured Low-Rank Matrix Completion"
% 2013, submitted to MRM. 
%
% It is recommended that this method be used to complete a calibration
% region and then use ESPIRiT to generate ESPIRiT maps. 
%
% INPUTS:
%       DATA [nx X ny X nc]  - 2D multi-coil k-space data with missing
%                               entries, preferably random! missing entries
%                               should be EXACTLY = 0!
%       kSize [ 1 x 2]       - sliding window size
%       wnRank [ 1 x 1]       - the window-normalized rank to enforce 
%                               (# of coils = full rank, 2.0 would be typical 
%                               for 8 coils with window size [5 5]
%       nIter                - number of iTerations.
%       show                 - show intermediate reconstruction show=0 will
%                               skip. Show = 100 will plot in figure 100
%
%
% Outputs:
%       res -                - 2D multi coil k-space data where the missing
%                               data was filled. 
%
% (c) Michael Lustig 2012
%



% mask = abs(kspace)>0;

% res = kspace;
% kspace_cart = gpuArray(kspace_cart);
N.S         = gpuArray(N.S);
kspace      = gpuArray(kspace);
sz = size(kspace, 3);

% [nsample, nos, nslice, ~, ncoil] = size(kspace);
for islice = 1:sz
    im(:, :, islice, :) = NUFFT.NUFFT_adj(kspace(:, :, islice, :), N); 
end


[sx, sy, sz, nc]	= size(im);
keep    = 1:floor(wnRank * prod(kernel_size(1:2)));
mask    = circular_mask_3d(sx, sx/2);
mask    = imresize(mask, [sx, sy]);
mask    = permute(mask, [1, 3, 2]);
mask    = imresize(mask, [sx, sz]);
mask    = permute(mask, [1, 3, 2]);
% us_mask = logical(abs(kspace(1, :, :, 1)));

% kspace_cart = kspace_cart .* mask;

% kspace_cart = gpuArray(kspace_cart);
% N.S         = gpuArray(N.S);
% kspace      = gpuArray(kspace);

for n = 1:nIter
    n
    tic
    % reorder data to get Hankel structure. 
    kspace_cart = fftshift2(fft2(fftshift2(im)));
    for islice = 1:sz
%         tic
        tmp = im2row(squeeze(kspace_cart(:, :, islice, :)), kernel_size(1:2));
        [tsx, tsy, tsz] = size(tmp);
        A = single(reshape(tmp, tsx, tsy*tsz));
        
        % SVD thresholding
        [U, S, V] = svdecon(A);
        A = U(:, keep) * S(keep,keep) * V(:,keep)';
        
        % Enforce Hankel structure
        A = reshape(A, tsx, tsy, tsz);
        tmp = single(row2im(A, [sx, sy, nc], kernel_size));
        
        % enforce data consistency
        im(:, :, islice, :) = permute(fftshift2(ifft2(fftshift2(tmp))), [1, 2, 4, 3]);
%         toc
    end
    toc
    
    figure(1)
    imagesc(abs(sos(im(:, :, 120, :))))
    axis image
    drawnow
    
    tic
    im = permute(im, [1, 3, 2, 4]);
    kspace_cart = fftshift2(fft2(fftshift2(im)));
    for islice = 1:sy
        tmp = im2row(squeeze(kspace_cart(:, :, islice, :)), kernel_size([1, 3]));
        [tsx, tsy, tsz] = size(tmp);
        A = single(reshape(tmp, tsx, tsy*tsz));
        
        % SVD thresholding
        [U, S, V] = svdecon(A);
        A = U(:, keep) * S(keep,keep) * V(:,keep)';
        
        % Enforce Hankel structure
        A = reshape(A, tsx, tsy, tsz);
        tmp = single(row2im(A, [sx, sz, nc], kernel_size));
        
        % enforce data consistency
        im(:, :, islice, :) = permute(fftshift2(ifft2(fftshift2(tmp))), [1, 2, 4, 3]);
    end
    im = permute(im, [1, 3, 2, 4]);
    toc
    
    figure(1)
    imagesc(abs(sos(im(:, :, 120, :))))
    axis image
    drawnow
    
    tic
    im = permute(im, [3, 2, 1, 4]);
    kspace_cart = fftshift2(fft2(fftshift2(im)));
    for islice = 1:sx
        tmp = im2row(squeeze(kspace_cart(:, :, islice, :)), kernel_size([1, 3]));
        [tsx, tsy, tsz] = size(tmp);
        A = single(reshape(tmp, tsx, tsy*tsz));
        
        % SVD thresholding
        [U, S, V] = svdecon(A);
        A = U(:, keep) * S(keep,keep) * V(:,keep)';
        
        % Enforce Hankel structure
        A = reshape(A, tsx, tsy, tsz);
        tmp = single(row2im(A, [sz, sy, nc], kernel_size));
        
        % enforce data consistency
        im(:, :, islice, :) = permute(fftshift2(ifft2(fftshift2(tmp))), [1, 2, 4, 3]);
        
    end
    im = permute(im, [3, 2, 1, 4]);
    toc
    
    tic
    for i = 1:sz
        f_update(:, :, i, :) = NUFFT.NUFFT(im(:, :, i, :), N);
    end
    f_update = permute(kspace, [1, 2, 4, 3]) - f_update;
    for i = 1:sz
        f_update_im(:, :, i, :) = NUFFT.NUFFT_adj(f_update(:, :, i, :), N);
    end
    
%     f_update = NUFFT.NUFFT_adj(kspace - NUFFT.NUFFT(im, N) .* us_mask, N);
    
    im = im + f_update_im * 0.5;
    clear f_update kspace_cart
%     kspace_cart = fftshift2(fft2(fftshift2(im))) .* mask;
    
%     res = tmp.*(1-mask) + kspace;
    figure(1)
    imagesc(abs(sos(im(:, :, 120, :))))
    axis image
    drawnow
    toc
end

%{
for n = 1:nIter
    tic
    % reorder data to get Hankel structure. 
    tmp = im2row(squeeze(kspace_cart), kernel_size); 
    [tsx, tsy, tsz] = size(tmp);
    A = reshape(tmp, tsx, tsy*tsz);
    
    % SVD thresholding
    [U, S, V] = svd(A, 'econ');
    
    A = U(:, keep) * S(keep,keep) * V(:,keep)';
    
    % Enforce Hankel structure
    A = reshape(A, tsx, tsy, tsz);
    tmp = row2im(A, [sx, sy, nc], kernel_size);
    
    % enforce data consistency
    s_update = permute(fftshift2(ifft2(fftshift2(tmp))), [1, 2, 4, 3]) - im;
    
    f_update = NUFFT.NUFFT_adj(kspace - NUFFT.NUFFT(im, N), N);
    
    im = im + f_update * 0.3 + s_update;
    kspace_cart = fftshift2(fft2(fftshift2(im))) .* mask;
    
%     res = tmp.*(1-mask) + kspace;
    figure(1)
    imagesc(sos(im))
    axis image
    drawnow
    toc
end

%}

