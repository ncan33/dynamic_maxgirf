function imo = TPS_yt(im,para)
% Im_out = TPS(Im_in,para)
% para.coor_source source points
% para.coor_target target points
% para.kernel 't' for thin plate and 'g' for gaussian

im_source = im;
coor_source_x = para.coor_source(:,2);
coor_source_y = para.coor_source(:,1);
coor_target_x = para.coor_target(:,2);
coor_target_y = para.coor_target(:,1);

kernel = para.kernel;

if kernel == 'g'
    r = para.kernel_width;
end

nol = size(coor_source_x,1);

[sx, sy] = size(im_source);

dx = bsxfun(@minus,coor_source_x',coor_source_x);
dy = bsxfun(@minus,coor_source_y',coor_source_y);
dr = sqrt(dx.^2+dy.^2);

if kernel == 'g'
    phi = exp(-dr.^2/r.^2);
elseif kernel == 't'
    phi = dr.^2.*log(dr); 
    phi(isnan(phi)) = 0;
end

temp_lu = [coor_source_x'; coor_source_y'; ones(1,nol)];
temp_rd = [coor_source_y, coor_source_x, ones(nol,1)];
B = [temp_lu zeros(3,3); phi temp_rd];
m_left = [B zeros(size(B)); zeros(size(B)) B];
m_right = [0; 0; 0; coor_target_x; 0; 0; 0; coor_target_y];
warning off
coef = m_left\m_right;

kx = coef(1:nol+3);
ky = coef(nol+4:end);

x = 1:sx;
y = 1:sy;

dx_all = bsxfun(@minus,x,coor_source_x);
dy_all = bsxfun(@minus,y,coor_source_y);
dy_all = permute(dy_all,[1 3 2]);
dr_all = sqrt(bsxfun(@plus,dx_all.^2,dy_all.^2));

if kernel == 'g'
    phi_all = exp(-dr_all.^2/r.^2);
elseif kernel == 't'
    phi_all = dr_all.^2.*log(dr_all);
    phi_all(isnan(phi_all)) = 0;
end

[X,Y] = meshgrid(1:sx,1:sy);
B = vertcat(phi_all,permute(Y,[3 2 1]),permute(X,[3 2 1]),ones([1 sx sy]));

x_out = squeeze(sum(bsxfun(@times,kx,B),1));
y_out = squeeze(sum(bsxfun(@times,ky,B),1));

imo = griddata(x_out,y_out,double(im_source),X',Y');

