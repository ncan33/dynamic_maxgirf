function ttv = apply_patch_cluster_ttv_3D(Image, patch)

n_patch = size(patch.idx);
n_cluster = size(patch.cluster, 2);
ttv = zeros(size(Image), 'like', Image);

for ipatch = 1:n_patch
    patch_temp = Image(patch.idx{ipatch});
    for icluster = 1:n_cluster
        if length(patch.cluster{ipatch, icluster}) > 2
%             ttv(patch.idx{ipatch}(:, :, :, patch.cluster{ipatch, icluster})) = ttv(patch.idx{ipatch}(:, :, :, patch.cluster{ipatch, icluster})) + patch_temp(:, :, :, patch.cluster{ipatch, icluster});
            ttv(patch.idx{ipatch}(:, :, :, patch.cluster{ipatch, icluster})) = ttv(patch.idx{ipatch}(:, :, :, patch.cluster{ipatch, icluster})) + compute_3DtTV_circ(patch_temp(:, :, :, patch.cluster{ipatch, icluster}), 1, eps('single'));
        end
    end
end

ttv = ttv .* patch.mask;