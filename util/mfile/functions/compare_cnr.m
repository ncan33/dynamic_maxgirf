function compare_cnr(img1, img2)

t = size(img1, 3);

temp = figure;
set(temp,'Position',[10,10,1000,1000]);
imagesc(abs(sum(img1,3)))
colormap gray
brighten(0.4)
axis image
axis off


for i=1:2
    h = imfreehand(gca,'closed',false); 
    mask = createMask(h);
    nop = sum(sum(mask));
    value1(:,i) = squeeze(sum(sum(mask.*abs(img1),1),2))/nop;
    value2(:,i) = squeeze(sum(sum(mask.*abs(img2),1),2))/nop;
end

Colors = colors;
figure
hold on

plot(value1(:,1), 'LineWidth', 2, 'LineStyle', '-', 'Color', Colors(1,:))
plot(value1(:,2), 'LineWidth', 2, 'LineStyle', '--', 'Color', Colors(1,:))

plot(value2(:,1), 'LineWidth', 2, 'LineStyle', '-','Color', Colors(2,:))
plot(value2(:,2), 'LineWidth', 2, 'LineStyle', '--','Color', Colors(2,:))

max_y = max([value1(:); value2(:)]);
min_y = min([value1(:); value2(:)]);
dy = max_y - min_y;
max_y = max_y + dy/10;
min_y = min_y - dy/10;

axis([1, t, min_y, max_y])

cnr1 = (value1(:,1)-value1(:,2))./value1(:,2);
cnr2 = (value2(:,1)-value2(:,2))./value2(:,2);

cnr1 = (value1(:,1)./value1(:,2));
cnr2 = (value2(:,1)./value2(:,2));

figure
hold on
plot(cnr1, 'LineWidth', 2);
plot(cnr2, 'LineWidth', 2);
title('CNR curve', 'FontSize', 26)

max_y = max([cnr1(:); cnr2(:)]);
min_y = min([cnr1(:); cnr2(:)]);
dy = max_y - min_y;
max_y = max_y + dy/10;
min_y = min_y - dy/10;

axis([1, t, min_y, max_y])

fprintf(sprintf('cnr 1 = %.2f\ncnr 2 = %.2f\n', mean(cnr1), mean(cnr2)))