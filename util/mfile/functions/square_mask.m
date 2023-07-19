function mask = square_mask(imsize, masksize)

if numel(imsize) == 1
    imsize = [imsize, imsize];
end

if numel(masksize) == 1
    masksize = [masksize, masksize];
end

mask = false(imsize);
mask(1:masksize(1), 1:masksize(2)) = true;

mask = circshift(mask, [floor((imsize(1) - masksize(1))/2+1), floor((imsize(2) - masksize(2))/2+1)]);