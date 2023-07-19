function kSpace_cart = pre_interp_radial_to_cart_oversample(kSpace_radial,x_coor,y_coor,over_sampling)

[sx,nor,nf,nc,nSMS,ns] = size(kSpace_radial);
sx_over = sx*over_sampling;
data = double(reshape(kSpace_radial,[sx,nor,nf*nc*nSMS*ns]));

[XI,YI] = meshgrid(-sx/2+1:sx/((sx+1)*over_sampling):sx/2);

kSpace_cart = single(zeros(round(sx_over ),round(sx_over ),nf*nc*nSMS*ns));

parfor k=1:nf*nc*nSMS*ns
    n = mod(k,nf);
    if n==0
        n = nf; 
    end
    
    X  = double(x_coor(:,:,n));
    Y  = double(y_coor(:,:,n));
    
    warning off
    
    kSpace_cart(:,:,k) = griddata(X,Y,data(:,:,k),XI,YI);
end

kSpace_cart(isnan(kSpace_cart)) = 0;
kSpace_cart = reshape(kSpace_cart, [sx_over,sx_over,nf,nc,nSMS,ns]);
%{
img = fftshift(ifft2(kSpace_cart));
xl = (sx_over -sx)/2+1;
xr = xl + sx -1;
img_cut = ifftshift(img(xl:xr,xl:xr,:));
kSpace_cut = fft2(img_cut);
kSpace_cut = reshape(kSpace_cut, [sx,sx,nf,nc,nSMS,ns]);
%}


end