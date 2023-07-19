function im = TFD(im)

im = circshift(im, [0, 0, -1]) - im;