function kSpace_cart = pre_interp_radial_to_cart_over(kSpace_radial,x_coor,y_coor,factor)

[sx,nor,nf,nc,nSMS,ns] = size(kSpace_radial);
sx_over = sx*factor;
data = double(reshape(kSpace_radial,[sx*nor,nf*nc*nSMS*ns]));

x = reshape(x_coor*factor,[sx*nor,nf]);
y = reshape(y_coor*factor,[sx*nor,nf]);

x = repmat(x,[1 nc*nSMS*ns]);
y = repmat(y,[1 nc*nSMS*ns]);

Xr = round(x);
Yr = round(y);
%[Xr,Yr] = meshgrid(-sx/2+1:sx/((sx+1)*factor):sx/2);

kSpace_cart = single(zeros(round(sx_over+1)*round(sx_over+1),nf*nc*nSMS*ns));
%kSpace_cart = single(zeros((sx+1)*(sx+1),nf*nc*nSMS*ns));
kSpace_r = single(zeros(sx*nor,nf*nc*nSMS*ns));

parfor k=1:nf*nc*nSMS*ns
    warning off
    kSpace_r(:,k) = griddata(x(:,k),y(:,k),data(:,k),Xr(:,k),Yr(:,k));
end

kSpace_r(isnan(kSpace_r)) = 0;

indx = sub2ind([sx_over+1,sx_over+1,nf*nc*nSMS*ns],Xr+sx_over/2+1,Yr+sx_over/2+1);

for i=1:nf*nc*nSMS*ns
    kSpace_cart(indx(:,i),i) = kSpace_r(:,i);
end

kSpace_cart = reshape(kSpace_cart,[sx_over+1,sx_over+1,nf,nc,nSMS,ns]);
kSpace_cart(1,:,:,:,:,:) = [];
kSpace_cart(:,1,:,:,:,:) = [];

img = fftshift(ifft2(kSpace_cart));
xl = (sx*factor-sx)/2+1;
xr = xl + sx -1;
img_cut = ifftshift(img(xl:xr,xl:xr,:));
kSpace_cut = fft2(img_cut);
kSpace_cut = reshape(kSpace_cut, [sx,sx,nf,nc,nSMS,ns]);

end