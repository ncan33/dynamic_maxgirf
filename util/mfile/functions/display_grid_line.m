function display_grid_line(displacement, Image_recon, Image_mask, nzp)

%% display initial grid lines
[sx, sy, nof] = size(Image_recon);
figure

grid_line_interval = 0.2 * nzp;
line_spacing = 5 * nzp;

grid_lines_horizontal_x = grid_line_interval:grid_line_interval:sx;
grid_lines_horizontal_y = 1.5*nzp : line_spacing : sy;
[grid_lines_horizontal_x, grid_lines_horizontal_y] = meshgrid(grid_lines_horizontal_x, grid_lines_horizontal_y);

grid_lines_vertical_x = 0.5*nzp : line_spacing : sx;
grid_lines_vertical_y = grid_line_interval:grid_line_interval:sy;
[grid_lines_vertical_x, grid_lines_vertical_y] = meshgrid(grid_lines_vertical_x, grid_lines_vertical_y);

% mask
[x, y] = meshgrid(1:sx, 1:sy);
mask_1 = Image_mask(:, :, 1);
horizontal_mask = interp2(x, y, single(mask_1), grid_lines_horizontal_x, grid_lines_horizontal_y, 'nearest');
vertical_mask = interp2(x, y, single(mask_1), grid_lines_vertical_x, grid_lines_vertical_y, 'nearest');
horizontal_mask(isnan(horizontal_mask)) = 0;
vertical_mask(isnan(vertical_mask)) = 0;

grid_lines_horizontal_x = grid_lines_horizontal_x(logical(horizontal_mask));
grid_lines_horizontal_y = grid_lines_horizontal_y(logical(horizontal_mask));

grid_lines_vertical_x = grid_lines_vertical_x(logical(vertical_mask));
grid_lines_vertical_y = grid_lines_vertical_y(logical(vertical_mask));

% draw
imagesc(Image_recon(:, :, 1))
axis image
axis off
colormap gray

hold on
plot(grid_lines_horizontal_x, grid_lines_horizontal_y, '.', 'Color', colors(1))
plot(grid_lines_vertical_x, grid_lines_vertical_y, '.', 'Color', colors(2))

%% deform grid lines
dx = squeeze(displacement(:, :, 1, :));
dy = squeeze(displacement(:, :, 2, :));

grid_line_horizontal_dx_dynamic = grid_lines_horizontal_x;
grid_line_horizontal_dy_dynamic = grid_lines_horizontal_y;
grid_line_vertical_dx_dynamic = grid_lines_vertical_x;
grid_line_vertical_dy_dynamic = grid_lines_vertical_y;


for i = 1:nof-1
    grid_line_horizontal_dx = interp2(x, y, dy(:, :, i), grid_line_horizontal_dx_dynamic(:, :, i), grid_line_horizontal_dy_dynamic(:, :, i));
    grid_line_horizontal_dx(isnan(grid_line_horizontal_dx)) = 0;
    grid_line_horizontal_dy = interp2(x, y, dx(:, :, i), grid_line_horizontal_dx_dynamic(:, :, i), grid_line_horizontal_dy_dynamic(:, :, i));
    grid_line_horizontal_dy(isnan(grid_line_horizontal_dy)) = 0;
    
    grid_line_vertical_dx = interp2(x, y, dy(:, :, i), grid_line_vertical_dx_dynamic(:, :, i), grid_line_vertical_dy_dynamic(:, :, i));
    grid_line_vertical_dx(isnan(grid_line_vertical_dx)) = 0;
    grid_line_vertical_dy = interp2(x, y, dx(:, :, i), grid_line_vertical_dx_dynamic(:, :, i), grid_line_vertical_dy_dynamic(:, :, i));
    grid_line_vertical_dy(isnan(grid_line_vertical_dy)) = 0;
    
    grid_line_horizontal_dx_dynamic(:, :, i + 1) = grid_line_horizontal_dx_dynamic(:, :, i) + grid_line_horizontal_dx;
    grid_line_horizontal_dy_dynamic(:, :, i + 1) = grid_line_horizontal_dy_dynamic(:, :, i) + grid_line_horizontal_dy;
    grid_line_vertical_dx_dynamic(:, :, i + 1) = grid_line_vertical_dx_dynamic(:, :, i) + grid_line_vertical_dx;
    grid_line_vertical_dy_dynamic(:, :, i + 1) = grid_line_vertical_dy_dynamic(:, :, i) + grid_line_vertical_dy;
    
    clf
    imagesc(Image_recon(:, :, i + 1)), axis image, colormap gray, axis off
    hold on
    plot(grid_line_horizontal_dx_dynamic(:, :, i + 1), grid_line_horizontal_dy_dynamic(:, :, i + 1), '.', 'Color', colors(1))
    plot(grid_line_vertical_dx_dynamic(:, :, i + 1), grid_line_vertical_dy_dynamic(:, :, i + 1), '.', 'Color', colors(2))
    title(i)
    drawnow
    pause
    
end
