function [csm] = csm_ss(img, aliasing, smoothing)
%
%   [csm] = ismrm_estimate_csm_walsh(img)
%
%   Estimates relative coil sensitivity maps from a set of coil images
%   using the eigenvector method described by Walsh et al. (Magn Reson Med
%   2000;43:682-90.)
%
%   INPUT:
%     - img     [x,y,coil]   : Coil images
%     - smooth  scalar       : Smoothing block size (defaults to 5)
%
%   OUTPUT:
%
%     - csm     [x,y,coil    : Relative coil sensitivity maps
%
%
%   Code is based on an original implementation by Peter Kellman, NHLBI,
%   NIH
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%

% The Walsh method is the same with SOS combination if not smoothing
% Ye 04/29/20
if nargin < 2
    smoothing = 5;
end

[sx, sy, ncoils] = size(img);

% normalize by root sum of squares magnitude
% mag = sqrt(sum(img .* conj(img),3))
mag = sos(img);
%s_raw = img ./ repmat(mag + eps,[1 1 ncoils]); clear mag; 
s_raw = bsxfun(@rdivide,img,mag+eps); clear mag;


% compute sample correlation estimates at each pixel location
% Rs = ismrm_correlation_matrix(s_raw);
Rs = permute(conj(s_raw), [1,2,4,3]).*(s_raw);
% compute correlation matrix for slice aliasing
aliasing = sum(squeeze(aliasing(:,:,:,1)), 4);
aliasing = aliasing./sos(aliasing);
Rsa = permute(conj(aliasing), [1,2,4,3]).*(aliasing);

% apply spatial smoothing to sample correlation estimates (NxN convolution)
if smoothing>1
	h_smooth = ones(smoothing)/(smoothing^2); % uniform smoothing kernel
    for m = 1:ncoils
        for n = 1:ncoils
            Rs(:,:,m,n) = conv2(Rs(:,:,m,n),h_smooth,'same');
            Rsa(:,:,m,n) = conv2(Rsa(:,:,m,n),h_smooth,'same');
        end
    end
end

for i=1:sx
    for j=1:sy
        Rs(i,j,:,:) = squeeze(Rsa(i,j,:,:))\squeeze(Rs(i,j,:,:));
    end
end

% compute dominant eigenvectors of sample correlation matrices
[csm,~] = ismrm_eig_power(Rs); % using power method

return 


%Utility functions provided by Peter Kellman, NIH.
function [Rs] = ismrm_correlation_matrix(s)
% function [Rs]=correlation_matrix(s);
%
% function correlation_matrix calculates the sample correlation matrix (Rs) for
% each pixel of a multi-coil image s(y,x,coil)
%
% input:
%    s   complex multi-coil image s(y,x,coil)
% output:
%    Rs  complex sample correlation matrices, Rs(y,x,coil,coil)

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

[rows,cols,ncoils] = size(s);
Rs = zeros(rows,cols,ncoils,ncoils,class(s)); % initialize sample correlation matrix to zero
for i=1:ncoils
    for j=1:i-1
		Rs(:,:,i,j)=s(:,:,i).*conj(s(:,:,j));
		Rs(:,:,j,i)=conj(Rs(:,:,i,j)); % using conjugate symmetry of Rs
    end
	Rs(:,:,i,i)=s(:,:,i).*conj(s(:,:,i));
end

return

function [v,d] = ismrm_eig_power(R)
% function [v,d]=eig_power(R);
%
% vectorized method for calculating the dominant eigenvector based on
% power method. Input, R, is an image of sample correlation matrices
% where: R(y,x,:,:) are sample correlation matrices (ncoil x ncoil) for each pixel
%
% v is the dominant eigenvector
% d is the dominant (maximum) eigenvalue

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

%rows=size(R,1);cols=size(R,2);ncoils=size(R,3);
[rows,cols,ncoils,~] = size(R);
N_iterations=2;
v=ones(rows,cols,ncoils); % initialize e.v.

d=zeros(rows,cols);
for i=1:N_iterations
    %v=squeeze(sum(R.*repmat(v,[1 1 1 ncoils]),3));
	v = squeeze(sum(bsxfun(@times,R,v),3));
%     d=ismrm_rss(v);
    d = sos(v);
    d( d < eps) = eps;
	%v=v./repmat(d,[1 1 ncoils]);
    v = bsxfun(@rdivide,v,d);
end

p1=angle(conj(v(:,:,1)));
% (optionally) normalize output to coil 1 phase
%v=v.*repmat(exp(sqrt(-1)*p1),[1 1 ncoils]);
v = bsxfun(@times,v,exp(1i*p1));
v = conj(v);

return

function y = ismrm_rss(x,dim)
%
%   [mag] = ismrm_rss(samples, dim)
%
%   Computes root-sum-of-squares along a single dimension.
%
%
%   INPUT:
%     - x   : multi-dimensional array of samples
%     - dim : dimension of reduction; defaults to last dimension
%
%   OUTPUT:
%
%     - y       : root sum of squares result
%
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip J. Beatty (philip.beatty@sri.utoronto.ca)
%


if nargin==1
    dim=ndims(x);
else
    if isempty(dim); dim=ndims(x); end
end

y = squeeze(sqrt(sum(real(x).^2 + imag(x).^2,dim)));
return

