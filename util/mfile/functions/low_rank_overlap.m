function Image_lr = low_rank_overlap(Image, patch_size, shift, tau)

[sx, sy, nof] = size(Image);


x_begin = 1:shift(1):sx;
y_begin = 1:shift(2):sy;

x_begin(x_begin + patch_size(1) - 1 > sx) = [];
y_begin(y_begin + patch_size(2) - 1 > sy) = [];

nx = length(x_begin);
ny = length(y_begin);
npatch = nx*ny;

[x_begin, y_begin] = meshgrid(x_begin, y_begin);
x_begin = x_begin(:);
y_begin = y_begin(:);

mask = zeros(size(Image), 'like', Image);
Image_lr = zeros(size(Image), 'like', Image);
for ipatch = 1:npatch
    idx_x = x_begin(ipatch) : x_begin(ipatch) + patch_size(1) - 1;
    idx_y = y_begin(ipatch) : y_begin(ipatch) + patch_size(2) - 1;
    patch_temp = Image(idx_x, idx_y, :);
    patch_temp = reshape(patch_temp, [prod(patch_size), nof]);
    [u, s, v] = svd(patch_temp, 0);
    s = s - tau;
    s(s<0) = 0;
    patch_temp = u * s * v';
    [idx_x, idx_y, idx_t] = meshgrid(idx_y, idx_x, 1:nof);
    patch_idx = sub2ind([sx, sy, nof], idx_y, idx_x, idx_t);
    mask(patch_idx) = mask(patch_idx) + 1;
    Image_lr(patch_idx) = Image_lr(patch_idx) + reshape(patch_temp, [patch_size, nof]);
end
mask(mask==0) = 1;
Image_lr = Image_lr ./ mask;
