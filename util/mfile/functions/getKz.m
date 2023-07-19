function kz_index = getKz

% generate weighted GA sampling order for 3dsos
% input:
%       N: number of inplane interleaves
%       M: number of planes
%       totlenth: number of interleaves
%       w: weighting factor of sampling
% ouput:
%       sked: sampling order
%
% This should be easily merged to psd.
% ZYH, 09/21/2014

M = 12;
totlength = 12*144;
w = 1;
GA = -(3-sqrt(5))/2; % reverse to use large GA

% weight = (1./[M/2-0.5:-1:0.5 0.5:1:M/2-0.5]).^w;
weight = (1./abs([1:M]-(M/2+0.5))).^w;
kz_index = zeros(1,totlength);

GAweight = cumsum(weight/sum(weight));

angles = angle(exp(1i*2*pi*(0:GA:GA*(totlength-1))))/2/pi+0.5;

% sort to kz
for n_M = M:-1:1
    kz_index(angles<=GAweight(n_M)) = n_M;
end
