function im = shrink_soft(im, weight)

im = sign(im) .* max(abs(im) - weight, 0);
