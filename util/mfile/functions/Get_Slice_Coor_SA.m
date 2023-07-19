function [Coor] = Get_Slice_Coor_SA(Isize,SliceInfo)

sx = Isize(1);
sy = Isize(2);
nsl = length(SliceInfo.SliceThickness);
[x,y] = meshgrid(1:sx,1:sy);

x = x-(sx+1)/2; 
y = y-(sx+1)/2;
Z = zeros(sx,sy);

Coor = zeros(sx*sy,nsl,3);

for i=1:nsl
    PixelSize_x = SliceInfo.FOV(1)/sx;
    PixelSize_y = SliceInfo.FOV(2)/sy;
    
    X = x * PixelSize_x;
    Y = y * PixelSize_y;
    coor = cat(2,X(:),Y(:),Z(:));
    
    dx = SliceInfo.SlicePosition(i,1);
    dy = SliceInfo.SlicePosition(i,2);
    dz = SliceInfo.SlicePosition(i,3);
    
    linear_shift = [dx,dy,dz];
    
    ix = SliceInfo.SliceOrintation(i,1);
    iy = SliceInfo.SliceOrintation(i,2);
    iz = SliceInfo.SliceOrintation(i,3);
    
    % construct rotation matrix along x
    theta_x = atan(iy/iz);
    Rx = [1 0 0; 0 cos(theta_x) -sin(theta_x); 0 sin(theta_x) cos(theta_x)];
    
    % rotate norm vector along x
    ii = sum(bsxfun(@times,Rx,[ix,iy,iz]),2);
    
    % construct rotation matrix along y
    theta_y = -atan(ii(1)/ii(3));
    Ry = [cos(theta_y) 0 sin(theta_y); 0 1 0; -sin(theta_y) 0 cos(theta_y)];
    
    % construct rotation matrix along z
    theta_z = SliceInfo.InPlaneRot(i);
    Rz = [cos(theta_z),-sin(theta_z),0;sin(theta_z),cos(theta_z),0;0,0,1];
    
    % combine rotation matrix x,y,z
    R = pinv(Rx)*pinv(Ry)*Rz;
    coor = bsxfun(@times,R,permute(coor,[3,2,1]));
    coor = squeeze(sum(coor,2));
    coor = coor.';

    Coor(:,i,:) = coor + linear_shift;
    
end
Coor = reshape(Coor,[sx,sy,nsl,3]);