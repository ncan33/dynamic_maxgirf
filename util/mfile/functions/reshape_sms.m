function img = reshape_sms(img)

img = squeeze(img);
[sx, sy, nof, nsms, nslice] = size(img);
img = permute(img, [1, 5, 2, 4, 3]);
img = reshape(img, [sx * nslice, sy * nsms, nof]);