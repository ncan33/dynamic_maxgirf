function display_local_max(Image, threshold)

if nargin == 1
    threshold = -inf;
end
[sx, sy, nof] = size(Image);
figure
for i = 1:nof
    clf
    im_frame = Image(:, :, i);
    imagesc(im_frame)
    colormap gray
    axis image
    regional_max = imregionalmax(im_frame, 4);
    regional_max((regional_max .* im_frame) < threshold) = false;
    [maxy, maxx] = find(regional_max);
    hold on
    plot(maxx, maxy, '.', 'Color', colors(1), 'MarkerSize', 15)
    
    title(i)
    pause
end