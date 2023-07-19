function display_strain(displacement, Image_recon, mask, nzp)

%% display initial grid lines
[sx, sy, nof] = size(Image_recon);


grid_line_interval = 1 * nzp;
line_spacing = 5 * nzp;

grid_lines_horizontal_x = grid_line_interval:grid_line_interval:sx;
grid_lines_horizontal_y = grid_line_interval:grid_line_interval:sx;
[grid_lines_horizontal_x, grid_lines_horizontal_y] = meshgrid(grid_lines_horizontal_x, grid_lines_horizontal_y);

% % draw
% figure
% imagesc(Image_recon(:, :, 1))
% axis image
% axis off
% colormap gray
% 
% hold on
% plot(grid_lines_horizontal_x, grid_lines_horizontal_y, '.', 'Color', colors(1))

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
    

%{
figure
for i = 1:nof-1
    clf
    ax1 = axes;
    imagesc(vec(x), vec(y), Image_recon(:, :, i + 1)), axis image, colormap gray, axis off
    P = get(ax1,'Position');
    XLIM = get(ax1,'XLim');
    YLIM = get(ax1,'YLim');
    PA = get(ax1,'PlotBoxAspectRatio');
    
    x_temp = x_dynamic(:, ceil(1+dd/2):end-floor(dd/2), i) .* mask(:, ceil(1+dd/2):end-floor(dd/2));
    y_temp = 100-y_dynamic(:, ceil(1+dd/2):end-floor(dd/2), i) .* mask(:, ceil(1+dd/2):end-floor(dd/2));
    strain_x_temp = -strain_x(:, :, i) .* mask(:, ceil(1+dd/2):end-floor(dd/2));
    
    ax2 = axes;
%     scatter(vec(x_dynamic(:, ceil(1+dd/2):end-floor(dd/2), i) .* mask(:, ceil(1+dd/2):end-floor(dd/2))), vec(100-y_dynamic(:, ceil(1+dd/2):end-floor(dd/2), i) .* mask(:, ceil(1+dd/2):end-floor(dd/2))), [], vec(-strain_x(:, :, i) .* mask(:, ceil(1+dd/2):end-floor(dd/2))), 'filled', 'square');
    surf(x_temp, y_temp, strain_x_temp, 'LineStyle', 'none');
    view(0, 90)
    
    linkaxes([ax1,ax2])
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    
    colormap(ax2,'jet')
    caxis([-0.4,0.4])
    
    set(ax2,'Position',P,'XLim',XLIM,'YLim',YLIM,'PlotBoxAspectRatio',PA)
    alpha(0.5)
    
    title(i)
    drawnow
    pause(0.1)
end
%}

file_name = 'strain_x.mp4';
make_video = 1;
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
    
    pause(0.1)
    
    if make_video
        ax = gca;
        ax.Units = 'pixels';
        pos = ax.Position;
        ti = ax.TightInset;
        
        rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
        image = getframe(ax,rect);

        writeVideo(v, image)
    end
end
close(v)





file_name = 'strain_y.mp4';
make_video = 1;
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
    
    pause(0.1)
    
    if make_video
        ax = gca;
        ax.Units = 'pixels';
        pos = ax.Position;
        ti = ax.TightInset;
        
        rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
        image = getframe(ax,rect);

        writeVideo(v, image)
    end
end
close(v)

%{
figure
for i = 1:nof-1
    clf
    ax1 = axes;
    imagesc(vec(x), vec(y), Image_recon(:, :, i + 1)), axis image, colormap gray, axis off
    P = get(ax1,'Position');
    XLIM = get(ax1,'XLim');
    YLIM = get(ax1,'YLim');
    PA = get(ax1,'PlotBoxAspectRatio');
    
    ax2 = axes;
    scatter(vec(x_dynamic(ceil(1+dd/2):end-floor(dd/2), :, i) .* mask(ceil(1+dd/2):end-floor(dd/2), :)), vec(100-y_dynamic(ceil(1+dd/2):end-floor(dd/2), :, i) .* mask(ceil(1+dd/2):end-floor(dd/2), :)), [], vec(-strain_y(:, :, i) .* mask(ceil(1+dd/2):end-floor(dd/2), :)), 'filled', 'square');
    
    linkaxes([ax1,ax2])
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    
    colormap(ax2,'jet')
    caxis([-0.4,0.4])
    
    
    set(ax2,'Position',P,'XLim',XLIM,'YLim',YLIM,'PlotBoxAspectRatio',PA)
    
    
    title(i)
    drawnow
    pause(0.1)
end
%}
