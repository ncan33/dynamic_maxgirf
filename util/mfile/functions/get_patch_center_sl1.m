function center_all = get_patch_center_sl1(input, patch_size)

[sx, sy, nc] = size(input);

nx = sx - patch_size(1) + 1;
ny = sy - patch_size(2) + 1;

patch_begin_x = 1:nx;
patch_begin_y = 1:ny;

patch_begin_x = patch_begin_x + round(patch_size(1)/2) - 1;
patch_begin_y = patch_begin_y + round(patch_size(2)/2) - 1;

center_all = zeros([nc, nx*ny], 'like', input);

idx = 0;
for i=1:nx
    for j=1:ny
        idx = idx + 1;
        x_idx = patch_begin_x(i);
        y_idx = patch_begin_y(j);
        center_all(:,idx) = input(x_idx, y_idx, :);
    end
end
