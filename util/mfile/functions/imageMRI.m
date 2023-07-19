function f = imageMRI(imshow)

flag = length(size(imshow)) == 4;

if flag
    [sx_, sy_, n1, n2] = size(imshow);
    imshow = permute(imshow, [1, 3, 2, 4]);
    imshow = reshape(imshow, sx_ * n1, sy_ * n2);
else
    imshow = imshow(:, :);
end

[sx, sy] = size(imshow);

f = figure;
if isreal(imshow)
    imagesc(imshow)
else
    imagesc(abs(imshow))
end
axis image
axis off
colormap gray

set(gcf, 'Position', [100, 100, sy / sx * 400, 400]);

set(gca, 'Pos', [0, 0, 1, 1]);

if flag
    hold on
    for i = 1 : n1-1
        plot([1, sy], [sx_, sx_] * i + 0.5, 'Color', colors(1), 'LineWidth', 2)
    end
    for i = 1 : n2-1
        plot([sy_, sy_] * i + 0.5, [1, sx], 'Color', colors(1), 'LineWidth', 2)
    end
end