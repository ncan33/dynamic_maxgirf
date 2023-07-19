function patchs = get_patch_cluster_3D(Image, patch_size, shift, n_cluster)

if ~isreal(Image)
    Image = abs(Image);
end
[sx, sy, sz, nof] = size(Image);


x_begin = 1:shift(1):sx;
y_begin = 1:shift(2):sy;
z_begin = 1:shift(3):sz;

x_begin(x_begin + patch_size(1) - 1 > sx) = [];
y_begin(y_begin + patch_size(2) - 1 > sy) = [];
z_begin(z_begin + patch_size(3) - 1 > sz) = [];

nx = length(x_begin);
ny = length(y_begin);
nz = length(z_begin);
npatch = nx*ny*nz;

[x_begin, y_begin, z_begin] = meshgrid(x_begin, y_begin, z_begin);
x_begin = x_begin(:);
y_begin = y_begin(:);
z_begin = z_begin(:);

patch = cell(npatch, 1);
cluster = cell(npatch, n_cluster);
mask = zeros(size(Image), 'like', Image);
for ipatch = 1:npatch
    idx_x = x_begin(ipatch) : x_begin(ipatch) + patch_size(1) - 1;
    idx_y = y_begin(ipatch) : y_begin(ipatch) + patch_size(2) - 1;
    idx_z = z_begin(ipatch) : z_begin(ipatch) + patch_size(3) - 1;
    patch_temp = Image(idx_x, idx_y, idx_z, :);
    patch_temp = reshape(patch_temp, [prod(patch_size), nof]).';
    clusters = kmeans(patch_temp, n_cluster);
    for icluster = 1:n_cluster
        cluster{ipatch, icluster} = find(clusters == icluster);
    end
    [idx_x, idx_y, idx_z, idx_t] = ndgrid(idx_y, idx_x, idx_z, 1:nof);
    patch{ipatch} = sub2ind([sx, sy, sz, nof], idx_y, idx_x, idx_z, idx_t);
    mask(patch{ipatch}) = mask(patch{ipatch}) + 1;
end
mask(mask==0) = 1;
mask = 1./mask;

patchs.idx = patch;
patchs.cluster = cluster;
patchs.mask = mask;