function mask = circular_mask(imsize, radius)

if numel(imsize) == 1
    imsize = [imsize, imsize];
end

mask = zeros(imsize);
mask(floor(imsize(1)/2 + 1), floor(imsize(2)/2 + 1)) = 1;
mask = bwdist(mask);
mask(mask>radius) = 0;
mask(floor(imsize(1)/2 + 1), floor(imsize(2)/2 + 1)) = 1;
mask = logical(mask);
