function display_optical_flow(displacement, image)

[sx, sy, nof] = size(image);
[x, y] = meshgrid(1:sx, 1:sy);

figure
interval = round(sx / 50);
for i = 1:nof-1
    clf
    imagesc(image(:,:,i))
    axis image
    axis off
    hold on
    quiver(x(20:interval:end-20, 20:interval:end-20), y(20:interval:end-20, 20:interval:end-20), displacement(20:interval:end-20, 20:interval:end-20,2,i), displacement(20:interval:end-20, 20:interval:end-20,1,i), 1.5, 'LineWidth', 2, 'Color', colors(3))
    colormap gray
    drawnow
    set(gca, 'pos', [0, 0, 1, 1])    
end