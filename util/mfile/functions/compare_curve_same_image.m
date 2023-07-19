function [value, mask_all] = compare_curve_same_image(Image,n)
if nargin == 1
    n = 1;
end

Image = (squeeze(Image));

temp = figure;
set(temp,'Position',[10,10,1000,1000]);
imagesc((sum(Image,3)))
colormap gray
brighten(0.4)
axis image
axis off

for i=1:n
    h = imfreehand(gca,'closed',false); 
    mask = createMask(h);
    nop = sum(sum(mask));
    value(:,i) = squeeze(sum(sum(mask.* Image,1),2))/nop;
    mask_all(:, :, i) = mask;
end

figure
hold on
for i=1:n
    plot(value(:,i))
    string(i,:) = ['curve',num2str(i)];
end
legend(string)
