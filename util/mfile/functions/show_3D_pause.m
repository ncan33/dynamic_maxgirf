function show_3D_pause(img, brightness, flag)
%--------------------------------------------------------------------------
%   show_3D_pause(img, brightness)
%--------------------------------------------------------------------------
%   Display 3D dynamic images frame-by-frame
%--------------------------------------------------------------------------
%   Inputs:      
%       - img           [sx, sy, sz, nof, ...]
%       - brightness    [scalar within 0-1]
%
%       - img           3D dynamic images with at least 4 dimensions
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%--------------------------------------------------------------------------

if nargin == 1
    brightness = 0.2;
    flag = 1;
elseif nargin == 2
    flag = 1;
end

if ~isreal(img)
    img = abs(img);
end
[sx,sy,ns,nf] = size(img);
ns_sqrt = ceil(sqrt(ns));
if ns_sqrt^2~=ns
    img(:,:,ns+1:ns_sqrt^2,:) = zeros(sx,sy,ns_sqrt^2-ns,nf);
    ns = ns_sqrt^2;
end
img = reshape(img,sx,sy*ns_sqrt,ns_sqrt,nf);
img = permute(img,[1 3 2 4]);
img = reshape(img,sx*ns_sqrt,sy*ns_sqrt,nf);
figure
for i=1:size(img,3)
    imagesc(img(:,:,i))
    colormap gray
    axis image
    brighten(brightness)
    title(i)
    drawnow
    if flag
        pause
    end
end
end