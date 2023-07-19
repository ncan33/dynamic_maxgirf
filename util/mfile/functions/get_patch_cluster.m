function patchs = get_patch_cluster(Image, patch_size, shift, n_cluster)

if ~isreal(Image)
    Image = abs(Image);
end
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

patch = cell(npatch, 1);
cluster = cell(npatch, n_cluster);
mask = zeros(size(Image), 'like', Image);

for ipatch = 1:npatch
    idx_x = x_begin(ipatch) : x_begin(ipatch) + patch_size(1) - 1;
    idx_y = y_begin(ipatch) : y_begin(ipatch) + patch_size(2) - 1;
    patch_temp = Image(idx_x, idx_y, :);
    patch_temp = reshape(patch_temp, [prod(patch_size), nof]).';
    % get kmeans seed
    [~, patch_temp_order] = sort(mean(patch_temp, 2));
    patch_temp_sort = patch_temp(patch_temp_order, :);
    init_ = round(1 : nof/n_cluster : nof);
    init_ = [init_, nof];
    for i = 1 : n_cluster
        initc(i, :) = mean(patch_temp_sort(init_(i):init_(i+1), :));
    end
    % kmeans clustering 
    clusters = kmeans(patch_temp, n_cluster, 'Start', initc);
    for icluster = 1:n_cluster
        cluster{ipatch, icluster} = find(clusters == icluster);
    end
    [idx_x, idx_y, idx_t] = meshgrid(idx_y, idx_x, 1:nof);
    patch{ipatch} = sub2ind([sx, sy, nof], idx_y, idx_x, idx_t);
    mask(patch{ipatch}) = mask(patch{ipatch}) + 1;
end

mask(mask==0) = 1;
mask = 1./mask;

patchs.idx = patch;
patchs.cluster = cluster;
patchs.mask = mask;