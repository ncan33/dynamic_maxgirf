function kernel = SSG_kernel_training_ye(source, target, options)
%--------------------------------------------------------------------------
%   kernel = SSG_kernel_training_ye(source, target, options)
%--------------------------------------------------------------------------
%   Train split slice GRAPPA or slice GRAPPA kernel
%--------------------------------------------------------------------------
%   Inputs:      
%       - source                    [sx, sy, nc, nsms]
%       - target                    [sx, sy, nc, nsms]
%       - options                   [structure]
%           options.patch_height    [positive integer]
%           options.patch_width     [positive integer]
%           options.alpha           [scalar]
%           options.SSG             [0 or 1]
%
%           'sx'    calibration k-space size in x dimension
%           'sy'    calibration k-space size in y dimension
%           'nc'    number of coils
%           'nsms'  number of SMS slices
%
%       - source                    input source k-space
%       - target                    calibration target k-space
%       - options                   training options
%           options.patch_height    kernel size y dimension
%           options.patch_width     kernel size x dimension
%           options.alpha           SSG weight
%           options.SSG             whether to use SSG (1) or SG (0)
%--------------------------------------------------------------------------
%   Output:
%       - kernel                    [h*w*nc, 1, nc, nsms]
%       
%           'h'     options.patch_height
%           'w'     options.patch_width
%           'nc'    number of coils
%           'nsms'  number of SMS slices
%
%       - kernel                    discription of output
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%--------------------------------------------------------------------------

h = options.patch_height;
w = options.patch_width;

alpha = options.alpha;

[ny,nx,nc,nz] = size(source);

% number of blocks
nb_x = nx-h+1;
nb_y = ny-w+1;

patch_s = zeros(h*w, nb_x*nb_y, nc, nz);
patch_t = zeros(h*w, nb_x*nb_y, nc, nz);

for c=1:nc
    for z=1:nz
        patch_s(:,:,c,z) = im2col(source(:, :, c, z), [h, w], 'sliding');
        patch_t(:,:,c,z) = im2col(target(:, :, c, z), [h, w], 'sliding');
    end
end

patch_s = reshape(patch_s, [h, w, nb_y, nb_x, nc, nz]);
patch_s = permute(patch_s,[3 4 6 1 2 5]);
S = reshape(patch_s,[nb_y*nb_x, nz, h*w*nc]);

patch_t = reshape(patch_t, [h, w, nb_y, nb_x, nc, nz]);
patch_t = permute(patch_t,[3 4 6 1 2 5]);
T = patch_t(:,:,:,ceil(h/2),ceil(w/2),:,:); % center point of each block
T = reshape(T,[nb_y*nb_x, nz, nc]);

if options.SSG ==1
    kernel = zeros([h*w*nc, 1, nc, nz]);
    
    for z = 1:nz
        S_SSG = zeros(size(S,3),size(S,3));
        for sl = 1:nz
            if z==sl
                S_SSG = S_SSG + (alpha^2)*squeeze(S(:,sl,:))'*squeeze(S(:,sl,:));
            else
                S_SSG = S_SSG + squeeze(S(:,sl,:))'*squeeze(S(:,sl,:));
            end
        end
        CMK = S_SSG\((alpha^2)*squeeze(S(:,z,:))');
        kernel(:,:,:,z) = CMK * squeeze(T(:,z,:)); 
    end
    
else
    S_SG = squeeze(sum(S,2));
    CMK = (S_SG'*S_SG)\S_SG';
    kernel = CMK * T(:,:);
    kernel = reshape(kernel, [h*w*nc, nz, nc]);
    kernel = permute(kernel, [1, 4, 3, 2]);
end

end



