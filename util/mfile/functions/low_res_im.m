function im_out = low_res_im(im,siz)

sx = size(im,1);
sy = size(im,2);
k = fft2(im);
k = fftshift2(k);
cut_x = (sx-siz(1))/2;
cut_y = (sy-siz(2))/2;
k = k(cut_x+1:cut_x+siz(1),cut_y+1:cut_y+siz(2),:,:,:,:,:);
k = fftshift2(k);
im_out = ifft2(k);