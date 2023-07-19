function [shading_map,sTV_mask] = get_coil_shading_map(Image)

[sx,sy,nof,ns,nsms] = size(Image);

Image = abs(Image(:,:,:));

shading_map = zeros(sx,sy,nof*ns*nsms,'like',Image);
sTV_mask = shading_map;

for i=1:nof*ns*nsms
    shading_map(:,:,i) = spline_2D_yt(Image(:,:,i));
    max_temp = max(max(crop_half_FOV(shading_map(:,:,i))));
    min_temp = min(min(crop_half_FOV(shading_map(:,:,i))));
    sTV_mask(:,:,i) = (shading_map(:,:,i)-min_temp)/(max_temp-min_temp)+1;
    sTV_mask(:,:,i) = (1./sTV_mask(:,:,i)-0.5)*2+1;
    sTV_mask(sTV_mask>2) = 2;
    shading_map(:,:,i) = shading_map(:,:,i)/max_temp;
end

sTV_mask = reshape(sTV_mask,[sx,sy,nof,ns,nsms]);
shading_map = reshape(shading_map,[sx,sy,nof,ns,nsms]);