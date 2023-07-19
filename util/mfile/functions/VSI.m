function VSI(Image)


[sx,sy,~] = size(Image);


kSpace = fft2(Image);

shift_x = 0.5;
shift_y = 0.5;
linear_phase_x = (-pi:2*pi/(sx-1):pi)*shift_x;
linear_phase_x = exp(1i*linear_phase_x);
linear_phase_x = fftshift(linear_phase_x);

linear_phase_y = (-pi:2*pi/(sx-1):pi)*shift_y;
linear_phase_y = exp(1i*linear_phase_y);
linear_phase_y = fftshift(linear_phase_y);

kSpace_1 = kSpace.*linear_phase_x;
kSpace_1 = kSpace_1.*linear_phase_y.';

Image_1 = ifft2(kSpace_1);


shift_x = -0.5;
shift_y = 0.5;
linear_phase_x = (-pi:2*pi/(sx-1):pi)*shift_x;
linear_phase_x = exp(1i*linear_phase_x);
linear_phase_x = fftshift(linear_phase_x);

linear_phase_y = (-pi:2*pi/(sx-1):pi)*shift_y;
linear_phase_y = exp(1i*linear_phase_y);
linear_phase_y = fftshift(linear_phase_y);

kSpace_1 = kSpace.*linear_phase_x;
kSpace_1 = kSpace_1.*linear_phase_y.';

Image_2 = ifft2(kSpace_1);


shift_x = 0.5;
shift_y = -0.5;
linear_phase_x = (-pi:2*pi/(sx-1):pi)*shift_x;
linear_phase_x = exp(1i*linear_phase_x);
linear_phase_x = fftshift(linear_phase_x);

linear_phase_y = (-pi:2*pi/(sx-1):pi)*shift_y;
linear_phase_y = exp(1i*linear_phase_y);
linear_phase_y = fftshift(linear_phase_y);

kSpace_1 = kSpace.*linear_phase_x;
kSpace_1 = kSpace_1.*linear_phase_y.';

Image_3 = ifft2(kSpace_1);


shift_x = -0.5;
shift_y = -0.5;
linear_phase_x = (-pi:2*pi/(sx-1):pi)*shift_x;
linear_phase_x = exp(1i*linear_phase_x);
linear_phase_x = fftshift(linear_phase_x);

linear_phase_y = (-pi:2*pi/(sx-1):pi)*shift_y;
linear_phase_y = exp(1i*linear_phase_y);
linear_phase_y = fftshift(linear_phase_y);

kSpace_1 = kSpace.*linear_phase_x;
kSpace_1 = kSpace_1.*linear_phase_y.';

Image_4 = ifft2(kSpace_1);

Image_VSI = zeros(2*sx,2*sy);

Image_VSI(1:2:end,1:2:end) = Image_4;
Image_VSI(2:2:end,1:2:end) = Image_2;
Image_VSI(1:2:end,2:2:end) = Image_3;
Image_VSI(2:2:end,2:2:end) = Image_1;

return
