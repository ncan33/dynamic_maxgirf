function [line_profile, line] = get_line_profile(im, line)


interp_factor = 4;


im = (squeeze(im));

if ~isreal(im)
    im = abs(im);
end

if nargin == 1
    temp = figure;
    set(temp,'Position',[10,10,800,800]);
    imagesc((sum(im,3)))
    colormap gray
    brighten(0.4)
    axis image
    axis off
    
    line = drawline;
    line = line.Position;
end

x1 = line(1, 1);
x2 = line(2, 1);
y1 = line(1, 2);
y2 = line(2, 2);

% length = abs(sqrt(x2.^2 + y2.^2) - sqrt(x1.^2 + y1.^2));
length = sqrt((x1 - x2).^2 + (y1 - y2).^2);
length = floor(length * interp_factor);

dx = (x2 - x1)/(length-1);
dy = (y2 - y1)/(length-1);

xx = x1:dx:x2;
yy = y1:dy:y2;


nframe = size(im, 3);
nx = size(im, 1);
ny = size(im, 2);

[x, y] = meshgrid(1:nx, 1:ny);

line_profile = zeros(length, nframe);

for i = 1:nframe
    line_profile(:, i) = interp2(x, y, im(:,:,i), xx, yy);
end

figure
imagesc(line_profile)
colormap gray
brighten(0.4)
% axis image

title 'Time-Line Profile'
xlabel 'Time'
ylabel 'Line'
set(gcf, 'position', [100, 100, 800, 200])
% set(gca, 'FontSize'
