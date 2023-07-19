function [m,mask] = get_mean_value(Image,n)
if nargin == 1
    n = 1;
end

%Image = abs(squeeze(Image));

figure
imagesc(Image)
colormap gray
brighten(0.4)
axis image
axis off

for i=1:n
    h = imfreehand(gca,'closed',false);
    mask(:,:,i)=createMask(h);
    nop = sum(sum(mask(:,:,i)));
    m(i) = squeeze(sum(sum(mask(:,:,i).*Image,1),2))/nop;
end
