function im = medfilt2_ye(im,filter_size)

siz = size(im);
if length(siz)~=2
    im = reshape(im,[siz(1:2),prod(siz(3:end))]);
else
    siz(3) = 1;
end

for i=1:prod(siz(3:end))
    im_real(:,:,i) = medfilt2(real(im(:,:,i)),filter_size);
    im_imag(:,:,i) = medfilt2(imag(im(:,:,i)),filter_size);
end
im = complex(im_real,im_imag);
im = reshape(im,siz);
end