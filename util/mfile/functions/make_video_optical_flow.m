function make_video_optical_flow(displacement, Image_recon, mask, nzp, file_name_pre)

make_video = 1;


[sx, sy, nof] = size(Image_recon);

grid_line_interval = 1 * nzp;
line_spacing = 5 * nzp;

grid_lines_horizontal_x = grid_line_interval:grid_line_interval:sx;
grid_lines_horizontal_y = grid_line_interval:grid_line_interval:sx;
[grid_lines_horizontal_x, grid_lines_horizontal_y] = meshgrid(grid_lines_horizontal_x, grid_lines_horizontal_y);

%% deform grid lines
dx = squeeze(displacement(:, :, 1, :));
dy = squeeze(displacement(:, :, 2, :));

x_dynamic = grid_lines_horizontal_x;
y_dynamic = grid_lines_horizontal_y;

[x, y] = meshgrid(1:sx, 1:sy);
for i = 1:nof-1
    dx_temp = interp2(x, y, dy(:, :, i), x_dynamic(:, :, i), y_dynamic(:, :, i));
    dx_temp(isnan(dx_temp)) = 0;
    dy_temp = interp2(x, y, dx(:, :, i), x_dynamic(:, :, i), y_dynamic(:, :, i));
    dy_temp(isnan(dy_temp)) = 0;
      
    x_dynamic(:, :, i + 1) = x_dynamic(:, :, i) + dx_temp;
    y_dynamic(:, :, i + 1) = y_dynamic(:, :, i) + dy_temp;
end

dd = 5;
l0 = nzp * dd;
strain_x = x_dynamic(:, 1+dd:end, :) - x_dynamic(:, 1:end-dd, :);
strain_x = (strain_x - l0) / l0;
strain_y = y_dynamic(1+dd:end, :, :) - y_dynamic(1:end-dd, :, :);
strain_y = (strain_y - l0) / l0;

for i = 1:nof
   d_strain_x(:, :, i) = griddata(x_dynamic(:, ceil(1+dd/2):end-floor(dd/2), i), y_dynamic(:, ceil(1+dd/2):end-floor(dd/2), i), strain_x(:, :, i), x, y);
   d_strain_y(:, :, i) = griddata(x_dynamic(ceil(1+dd/2):end-floor(dd/2), :, i), y_dynamic(ceil(1+dd/2):end-floor(dd/2), :, i), strain_y(:, :, i), x, y);
   d_mask(:, :, i) = griddata(x_dynamic(:, :, i), y_dynamic(:, :, i), double(mask), x, y);
end
d_strain_x(isnan(d_strain_x)) = 0;
d_strain_y(isnan(d_strain_y)) = 0;
d_mask(isnan(d_mask)) = 0;

%% horizontal strain
file_name = sprintf('%s_strain_x.mp4', file_name_pre);

v = VideoWriter(file_name, 'MPEG-4');
v.FrameRate = 1000/5.58/2; % frames per second
v.Quality = 100;
open(v)

figure

for i = 1:nof-1
    clf
    ax2 = axes;
    x_temp = x_dynamic(:, ceil(1+dd/2):end-floor(dd/2), i);
    y_temp = 100-y_dynamic(:, ceil(1+dd/2):end-floor(dd/2), i);
    x_temp(x_temp>sx) = sx;
    y_temp(y_temp>sy) = sy;
    strain_x_temp = strain_x(:, :, i); % .* mask(:, ceil(1+dd/2):end-floor(dd/2));
    surf(x_temp, y_temp, strain_x_temp, 'LineStyle', 'none');
    view(0, 90)
    colormap(ax2,'jet')
    caxis([-0.4, 0.4])
    colorbar
    set(gca, 'FontSize', 20)
    
    ax1 = axes;
    I = imagesc(vec(x), vec(y), Image_recon(:, :, i + 1)); axis image, colormap(ax1,'gray'), axis off
    P = get(ax1,'Position');
    XLIM = get(ax1,'XLim');
    YLIM = get(ax1,'YLim');
    PA = get(ax1,'PlotBoxAspectRatio');
    
    linkaxes([ax1,ax2])
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    
    set(ax2,'Position',P,'XLim',XLIM,'YLim',YLIM,'PlotBoxAspectRatio',PA)
    
    set(I, 'AlphaData', 1 - d_mask(:, :, i) * 0.4);
    title('Horizontal Strain', 'FontSize', 26)
    
    set(gcf, 'Position', [20, 200, 650, 470])
    drawnow

    if make_video
        ax = gca;
        ax.Units = 'pixels';
        pos = ax.Position;
        ti = ax.TightInset;
        
        rect = [-ti(1)-pos(1)/2, -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
        image = getframe(ax,rect);

        writeVideo(v, image)
    end
end
close(v)



%% vertical strain

file_name = sprintf('%s_strain_y.mp4', file_name_pre);
v = VideoWriter(file_name, 'MPEG-4');
v.FrameRate = 1000/5.58/2; % frames per second
v.Quality = 100;
open(v)

figure

for i = 1:nof-1
    clf
    ax2 = axes;
    x_temp = x_dynamic(ceil(1+dd/2):end-floor(dd/2), :, i);
    y_temp = 100 - y_dynamic(ceil(1+dd/2):end-floor(dd/2), :, i);
    x_temp(x_temp>sx) = sx;
    y_temp(y_temp>sy) = sy;
    strain_y_temp = strain_y(:, :, i); % .* mask(:, ceil(1+dd/2):end-floor(dd/2));
    surf(x_temp, y_temp, strain_y_temp, 'LineStyle', 'none');
    view(0, 90)
    colormap(ax2,'jet')
    caxis([-0.4, 0.4])
    colorbar
    set(gca, 'FontSize', 20)
    
    ax1 = axes;
    I = imagesc(vec(x), vec(y), Image_recon(:, :, i + 1)); axis image, colormap(ax1,'gray'), axis off
    P = get(ax1,'Position');
    XLIM = get(ax1,'XLim');
    YLIM = get(ax1,'YLim');
    PA = get(ax1,'PlotBoxAspectRatio');
    
    linkaxes([ax1,ax2])
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    
    set(ax2,'Position',P,'XLim',XLIM,'YLim',YLIM,'PlotBoxAspectRatio',PA)
    
    set(I, 'AlphaData', 1 - d_mask(:, :, i) * 0.4);
    title('Vertical Strain', 'FontSize', 26)
    
    set(gcf, 'Position', [20, 200, 650, 470])
    drawnow
    
    if make_video
        ax = gca;
        ax.Units = 'pixels';
        pos = ax.Position;
        ti = ax.TightInset;
        
        rect = [-ti(1)-pos(1)/2, -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
        image = getframe(ax,rect);

        writeVideo(v, image)
    end
end
close(v)

%% 


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
mask_1 = mask(:, :, 1);
horizontal_mask = interp2(x, y, single(mask_1), grid_lines_horizontal_x, grid_lines_horizontal_y, 'nearest');
vertical_mask = interp2(x, y, single(mask_1), grid_lines_vertical_x, grid_lines_vertical_y, 'nearest');
horizontal_mask(isnan(horizontal_mask)) = 0;
vertical_mask(isnan(vertical_mask)) = 0;

grid_lines_horizontal_x = grid_lines_horizontal_x(logical(horizontal_mask));
grid_lines_horizontal_y = grid_lines_horizontal_y(logical(horizontal_mask));

grid_lines_vertical_x = grid_lines_vertical_x(logical(vertical_mask));
grid_lines_vertical_y = grid_lines_vertical_y(logical(vertical_mask));


%% deform grid lines
dx = squeeze(displacement(:, :, 1, :));
dy = squeeze(displacement(:, :, 2, :));

grid_line_horizontal_dx_dynamic = grid_lines_horizontal_x;
grid_line_horizontal_dy_dynamic = grid_lines_horizontal_y;
grid_line_vertical_dx_dynamic = grid_lines_vertical_x;
grid_line_vertical_dy_dynamic = grid_lines_vertical_y;


file_name = sprintf('%s_grid_line.mp4', file_name_pre);
make_video = 1;
v = VideoWriter(file_name, 'MPEG-4');
v.FrameRate = 1000/5.58/2; % frames per second
v.Quality = 100;
open(v)

figure

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
    
    
    title('Tagging Lines Deformation', 'FontSize', 26)
    
    set(gcf, 'Position', [20, 200, 650, 470])
    drawnow
    
    if make_video
        ax = gca;
        ax.Units = 'pixels';
        pos = ax.Position;
        ti = ax.TightInset;
        rect = [-ti(1)-pos(1)/2, -ti(2), pos(3)+ti(1)+ti(3)-pos(1)/2, pos(4)+ti(2)+ti(4)];
        
        image = getframe(ax, rect);

        writeVideo(v, image)
    end
    
end
close(v)




%% quiver
dd = displacement;

[x, y] = meshgrid(1:sx, 1:sy);
dd = dd .* permute(d_mask, [1, 2, 4, 3]);
nof = size(dd, 4);

arrow_interval = 1;

xx = x(1:arrow_interval:end, 1:arrow_interval:end);
yy = y(1:arrow_interval:end, 1:arrow_interval:end);

dx = squeeze(dd(1:arrow_interval:end, 1:arrow_interval:end, 1, :));
dy = squeeze(dd(1:arrow_interval:end, 1:arrow_interval:end, 2, :));

mask_small = logical(d_mask(1:arrow_interval:end, 1:arrow_interval:end, :));

max_d_all = max(vec(sqrt(dx.^2 + dy.^2)));

cmap = colormap(hsv(360));

mask2 = ones(sx, sy);
mask2(1:2:end, :) = 0;
mask2(:, 1:2:end) = 0;
mask2 = logical(mask2);

file_name = sprintf('%s_quiver.mp4', file_name_pre);
v = VideoWriter(file_name,'MPEG-4');
v.FrameRate = 1000/5.58 /2; % frames per second
v.Quality = 100;
open(v)

figure
for i = 1:nof-1
    clf
    imagesc(Image_recon(:,:,i))
    axis image
    axis off
    hold on
    colormap gray
    
    mask_temp = mask_small(:, :, i) & mask2;
    dx_temp = dx(:, :, i);
    dx_temp = dx_temp(mask_temp);
    
    dy_temp = dy(:, :, i);
    dy_temp = dy_temp(mask_temp);
    
    x_temp = xx(mask_temp);
    y_temp = yy(mask_temp);
    
    max_d = max(sqrt(dx_temp.^2 + dy_temp.^2));
%     dx_temp = dx_temp ./ max_d;
%     dy_temp = dy_temp ./ max_d;
    
    gloable_scale = max_d / max_d_all;
    
    theta = atan2(dx_temp, dy_temp);
    r = sqrt(dx_temp.^2 + dy_temp.^2);
%     r = (r + 0.1) / 1.1;
    
    dx_temp = r .* sin(theta);
    dy_temp = r .* cos(theta);
    
    for j = 1:length(dy_temp)
        theta = atan2(dx_temp(j), dy_temp(j));
        idx = ceil(rad2deg(theta));
        if theta <= 0 
            idx = idx + 360;
        end
        current_color = cmap(idx,:);
        relative_length = sqrt(dy_temp(j).^2 + dx_temp(j).^2) / max_d_all;
        alpha = (relative_length + 0.5) / 1.5;
        
        plot_arrow_(x_temp(j), y_temp(j), x_temp(j) + dy_temp(j) * 20, y_temp(j) + dx_temp(j) * 20, current_color, 1, 1.5);
    end
    drawnow
    title('Quiver Motion', 'FontSize', 26)
    set(gcf, 'Position', [20, 200, 650, 470])
    
    plot_circle([16, 85], 10, cmap)
    
    if make_video
        ax = gca;
        ax.Units = 'pixels';
        pos = ax.Position;
        ti = ax.TightInset;
        rect = [-ti(1)-pos(1)/2, -ti(2), pos(3)+ti(1)+ti(3)-pos(1)/2, pos(4)+ti(2)+ti(4)];
        
        image = getframe(ax, rect);

        writeVideo(v, image)
    end
end

close(v)