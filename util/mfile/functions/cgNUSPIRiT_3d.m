function [res,FLAG,RELRES,ITER,RESVEC,LSVEC] = cgNUSPIRiT_3d(kspace, x0, N, kernel, niter, lambda)

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
sx      = size(x0, 1);
sy      = size(x0, 2);
sz      = size(x0, 3);
x0      = x0 * sqrt(sx * sy * sz);

undersample_mask = logical(abs(kspace(1, :, :, :, 1)));

%%
n           = numel(x0);
imSize      = size(x0);
dataSize    = size(kspace);

b = [kspace(:); zeros([prod(imSize), 1])];

kernel  = double(kernel);
b       = double(b);
x0      = double(x0);

clear kspace 

[res, FLAG, RELRES, ITER, RESVEC, LSVEC] = lsqr(@(x,tflag)afun(x, N, kernel, dataSize, imSize, lambda, tflag, undersample_mask), b, [], niter, speye(n,n), speye(n,n), x0(:));

res = reshape(res, imSize);

function [y, tflag] = afun(x, N, kernel_image, dataSize, imSize, lambda, tflag, undersample_mask)

if strcmp(tflag, 'transp')
    tic
    fprintf('transpose \n')
    x1 = reshape(x(1:prod(dataSize)), dataSize);
    x2 = reshape(x(prod(dataSize)+1:end), imSize);
    
    y1 = zeros(imSize);
    for iz = 1:imSize(3)
        y1(:, :, iz, :) = NUFFT.NUFFT_adj(x1(:, :, iz, :), N);
    end
    
    y1 = fftshift(y1, 3);
    y1 = ifft(y1, [], 3);
    y1 = fftshift(y1, 3);
    
    y1 = y1 * sqrt(prod(imSize(1:3)));
    
    y2 = sum(conj(kernel_image) .* permute(x2, [1, 2, 3, 5, 4]), 5);
    y2 = y2 - x2;
    
    y = y1 + lambda * y2;
    y = y(:);
    toc
    
    %    x1 = double(NUFFT.NUFFT_adj(x1, N)) * sqrt(size(x1, 1) * size(x1, 2));
    %    x2_ = sum(conj(kernel_image) .* x2, 4);
    %    x2_ = permute(x2_, [1, 2, 4, 3]);
    %    x2_ = x2_ - x2;
    
    %    y = x1 + lambda * x2_;
    %    y = y(:);
    %    y = N'.* x1 + lambda * (kernel_image'*x2);
    %    y = y(:);
else
    tic
    
    fprintf('forward \n')
    x = reshape(x, imSize);
    
    y2 = sum(kernel_image .* x, 4);
    y2 = squeeze(y2) - x;
    
    x = fftshift(x, 3);
    x = fft(x, [], 3);
    x = fftshift(x, 3);
    y1 = zeros(dataSize);
    for iz = 1:dataSize(3)
        y1(:, :, iz, :, :) = NUFFT.NUFFT(x(:, :, iz, :), N);
    end
    y1 = y1 .* undersample_mask / sqrt(size(x, 1) * size(x, 2) * size(x, 3)) .* double(N.W);
    
    %     y1 = double(NUFFT.NUFFT(x / sqrt(size(x, 1) * size(x, 2) * size(x, 3)), N) .* undersample_mask .* N.W);
    %     y2 = sum(permute(kernel_image, [1, 2, 4, 3]) .* x, 4);
    %     y2 = y2 - squeeze(x);
    
    y = [y1(:); lambda * y2(:)];
    toc
end

