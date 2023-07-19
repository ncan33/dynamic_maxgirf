function Image_out = zero_pad_image(Image, inputsiz)

siz = size(Image);
if ~exist('inputsiz')
    inputsiz = siz(1:2)*2;
end

nzeros = round((inputsiz-siz(1:2))/2);
k = fftshift2(fft2(Image));
k = cat(1,zeros([nzeros(1),siz(2:end)],'like',Image),k,zeros([nzeros(1),siz(2:end)],'like',Image));
k = cat(2,zeros([inputsiz(1),nzeros(2),siz(3:end)],'like',Image),k,zeros([inputsiz(1),nzeros(2),siz(3:end)],'like',Image));
k = fftshift2(k);
Image_out = ifft2(k);
Image_out = Image_out*prod(inputsiz)/prod(siz(1:2));
if isreal(Image)
    Image_out = abs(Image_out);
end
