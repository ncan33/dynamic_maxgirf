function center_all = get_patch_center_section_sl1(input, patch_size)

[sx, sy, nc] = size(input);
mask = input == 0 ;

nx = sx - patch_size(1) + 1;
ny = sy - patch_size(2) + 1;

patch_begin_x = 1:nx;
patch_begin_y = 1:ny;

begin_x = patch_begin_x + round(patch_size(1)/2) - 1;
begin_y = patch_begin_y + round(patch_size(2)/2) - 1;

patch_end_x = patch_begin_x + patch_size(1) - 1;
patch_end_y = patch_begin_y + patch_size(2) - 1;

center_all = zeros([nc, 1], 'like', input);

idx = 0;
for i=1:nx
    for j=1:ny
        x_patch_idx = patch_begin_x(i):patch_end_x(i);
        y_patch_idx = patch_begin_y(j):patch_end_y(j);
        idx_temp = sum(vec(mask(x_patch_idx, y_patch_idx, 1)));
        if ~idx_temp
            idx = idx + 1;
            x_idx = begin_x(i);
            y_idx = begin_y(j);
            center_all(:,idx) = input(x_idx, y_idx, :);
        end
    end
end
