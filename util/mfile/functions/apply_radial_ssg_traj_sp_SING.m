function out = apply_radial_ssg_traj_sp_SING(kSpace_mb, idx, kernel, nsms)

[sx, nor, nof, nc] = size(kSpace_mb);

np = length(idx);
out = zeros([np, nc, nsms], 'like', kSpace_mb);

kSpace_mb = reshape(kSpace_mb, [sx*nor*nof, nc]);
for i=1:np
    source = kSpace_mb(idx{i}, :);
    target = kernel{i} * source(:);
    out(i, :, :) = reshape(target, [nc, nsms]);
end