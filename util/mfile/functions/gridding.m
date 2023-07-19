function image = gridding(kSpace, kx, ky, im_size)

[~, nor, nof, nc, nslice] = size(kSpace);
image = zeros([im_size, nof, nc, nslice], 'like', kSpace);

for ix = 1:im_size(1)
    ix
    for iy = 1:im_size(2)
        image(ix, iy, :, :, :) = sum(sum(kSpace .* exp(2*pi*1i .* (kx + ky))));
    end
end