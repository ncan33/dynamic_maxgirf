function Image_lr = apply_patch_cluster_low_rank(Image, patch, tau)

n_patch = size(patch.idx);
patch_x = size(patch.idx{1}, 1);
patch_y = size(patch.idx{1}, 2);
n_cluster = size(patch.cluster, 2);
Image_lr = zeros(size(Image), 'like', Image);

for ipatch = 1:n_patch
    patch_temp = Image(patch.idx{ipatch});
    for icluster = 1:n_cluster
        cluster_length = length(patch.cluster{ipatch, icluster});
        patch_temp_c = reshape(patch_temp(:,:,patch.cluster{ipatch, icluster}), [patch_x * patch_y, cluster_length]);
        [u, s, v] = svd(patch_temp_c, 0);
        s = s - tau;
        s(s<0) = 0;
        patch_temp_c = u * s * v';
        
        Image_lr(patch.idx{ipatch}(:, :, patch.cluster{ipatch, icluster})) = Image_lr(patch.idx{ipatch}(:, :, patch.cluster{ipatch, icluster})) + reshape(patch_temp_c, [patch_x, patch_y, cluster_length]);
    end
end

Image_lr = Image_lr .* patch.mask;