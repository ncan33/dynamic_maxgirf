function [curve] = compare_curve(Image1,Image2, n)
if nargin == 2
    n = 1;
end

Image1 = (squeeze(Image1));
Image2 = (squeeze(Image2));

Image1 = Image1(:,:,:,1);
Image2 = Image2(:,:,:,1);

figure
imagesc(sum(Image1,3))
colormap gray
brighten(0.4)
axis image
axis off

for i=1:n
    h = imfreehand(gca,'closed',false);
    mask=createMask(h);
    nop = sum(sum(mask));
    value1(:, i) = squeeze(sum(sum(mask.*Image1,1),2))/nop;
    value2(:, i) = squeeze(sum(sum(mask.*Image2,1),2))/nop;
end
%close


%value2 = value2/mean(value2(:))*mean(value1(:));
figure
hold on
if size(value1) == 1
    for i=1:n
        plot(i,value1(:, i),'*')
        plot(i,value2(:, i),'*')
    end
else
    for i=1:n
        plot(value1(:, i), 'Color', colors(i), 'LineWidth', 2)
        plot(value2(:, i), 'LineStyle', '--','Color', colors(i), 'LineWidth', 2)
    end
end
legend('curve 1', 'curve2')

curve = [value1,value2];