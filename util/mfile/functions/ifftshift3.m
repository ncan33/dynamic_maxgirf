function in = ifftshift3(in)
%--------------------------------------------------------------------------
%   [in] = ifftshift3(in)
%--------------------------------------------------------------------------
%   Perform 3D ifftshift
%--------------------------------------------------------------------------
%   Inputs:      
%       - in    [sx, sy, sz, ...]
%
%       - in    array with at least 3 dimensions
%--------------------------------------------------------------------------
%   Output:
%       - in    [sx, sy, sz, ...]
%
%       - in    3D ifftshift of input
%--------------------------------------------------------------------------
%   Author:
%       Ye Tian
%       E-mail: phye1988@gmail.com
%--------------------------------------------------------------------------

in = ifftshift(in,1);
in = ifftshift(in,2);
in = ifftshift(in,3);