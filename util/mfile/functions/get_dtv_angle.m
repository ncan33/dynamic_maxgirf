function phi = get_dtv_angle(Image)

phi = zeros(size(Image),'like',Image);
phi(1:end-1,:,:,:,:,:,:) = abs(diff(Image,1,1));
phi(:,1:end-1,:,:,:,:,:) = phi(:,1:end-1,:,:,:,:,:) + 1i*abs(diff(Image,1,2));
phi = angle(phi);