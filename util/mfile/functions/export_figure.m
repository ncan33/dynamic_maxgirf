function export_figure(image, imgname)

f = figure;
imagesc(image)
axis image
axis off
colormap gray
set(gcf, 'Position', [0, 0, 400, 400])
set(gca, 'Pos', [0 0 1 1])
set(gca, 'FontSize', 20)
% brighten(0.2)
hgexport(f, imgname)

