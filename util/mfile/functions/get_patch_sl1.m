function patch_all = get_patch_sl1(input, patch_size)

[sx, sy, nc] = size(input);

nx = sx - patch_size(1) + 1;
ny = sy - patch_size(2) + 1;

patch_begin_x = 1:nx;
patch_begin_y = 1:ny;

patch_end_x = patch_begin_x + patch_size(1) - 1;
patch_end_y = patch_begin_y + patch_size(2) - 1;

patch_all = zeros([patch_size, nc, nx*ny], 'like', input);

idx = 0;
for i=1:nx
    for j=1:ny
        idx = idx + 1;
        x_idx = patch_begin_x(i):patch_end_x(i);
        y_idx = patch_begin_y(j):patch_end_y(j);
        patch_all(:,:,:,idx) = input(x_idx, y_idx, :);
    end
end

patch_all = reshape(patch_all, [prod(patch_size)*nc, idx]);