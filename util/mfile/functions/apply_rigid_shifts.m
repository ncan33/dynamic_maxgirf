function im = apply_rigid_shifts(im, dx, dy)

[sx, sy] = size(im);

phasex = reshape(-2i*pi*((1:sx)-floor(sx/2)-1)/sx,sx,1);
phasey = reshape(-2i*pi*((1:sy)-floor(sy/2)-1)/sy,1,sy);
phase = exp(phasex*dx + phasey*dy);

im = fftshift2(ifft2(fftshift2(fftshift2(fft2(fftshift2(im))).*phase)));
