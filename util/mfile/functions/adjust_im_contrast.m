function im = adjust_im_contrast(im, pmin, pmax)

pmin_input = min(vec(im));
pmax_input = max(vec(im));

im = (im - pmin_input) / (pmax_input - pmin_input);

im(im<pmin) = pmin;
im(im>pmax) = pmax;
im = (im - pmin) / (pmax - pmin);