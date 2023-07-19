function [res,FLAG,RELRES,ITER,RESVEC,LSVEC] = cgNUSPIRiT_2d(kspace, x0, N, kernel_image, niter, lambda)

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
x0      = x0 * sqrt(size(x0, 1) * size(x0, 2));

undersample_mask = logical(abs(kspace(1, :, :, 1)));

%%
n           = numel(x0);
imSize      = size(x0);
dataSize    = size(kspace);

b = [kspace(:); zeros([prod(imSize), 1], 'single')];

kernel_image = double(kernel_image);
b = double(b);
x0 = double(x0);

tic
[res, FLAG, RELRES, ITER, RESVEC, LSVEC] = lsqr(@(x,tflag)afun(x, N, kernel_image, dataSize, imSize, lambda, tflag, undersample_mask), b, [], niter, speye(n,n), speye(n,n), x0(:));
toc

res = reshape(res, imSize);

function [y, tflag] = afun(x, N, kernel_image, dataSize, imSize, lambda, tflag, undersample_mask)

if strcmp(tflag, 'transp')
   x1 = reshape(x(1:prod(dataSize)), dataSize);
   x2 = reshape(x(prod(dataSize)+1:end), imSize);
   
   x1 = double(NUFFT.NUFFT_adj(x1, N)) * sqrt(size(x1, 1) * size(x1, 2));
   x2_ = sum(conj(kernel_image) .* x2, 4);
   x2_ = permute(x2_, [1, 2, 4, 3]);
   x2_ = x2_ - x2;
   
   y = x1 + lambda * x2_;
   y = y(:);
%    y = N'.* x1 + lambda * (kernel_image'*x2);
%    y = y(:);
else
    x = reshape(x, imSize);
    y1 = double(NUFFT.NUFFT(x / sqrt(size(x, 1) * size(x, 2)), N) .* undersample_mask .* N.W);
    y2 = sum(permute(kernel_image, [1, 2, 4, 3]) .* x, 4);
    y2 = y2 - squeeze(x);
%     y1 = N .* x;
%     y2 = kernel_image * x;
    y = [y1(:); lambda * y2(:)];
end

