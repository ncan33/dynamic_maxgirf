function dxy = get_total_variation(Image)

dx = diff(Image,1,1);
dy = diff(Image,1,2);

dx(end+1,:,:,:,:,:,:) = 0;
dy(:,end+1,:,:,:,:,:) = 0;

dxy = (dx.^2 + dy.^2).^0.5;