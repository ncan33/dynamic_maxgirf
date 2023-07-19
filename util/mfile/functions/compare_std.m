function compare_std(Image1,Image2)

Image1 = abs(squeeze(Image1));
Image2 = abs(squeeze(Image2));

nof = size(Image1,3);

Image1 = Image1(:,:,:,1);
Image2 = Image2(:,:,:,1);

figure
imagesc(sum(Image1,3))
colormap gray
brighten(0.4)
axis image
axis off

h = imfreehand(gca,'closed',false); 
mask=createMask(h);
nop = sum(sum(mask));
%close
value1 = mask.*Image1;
value1 = value1(value1~=0);
value1 = reshape(value1,length(value1)/nof,nof);
std1 = std(value1);
value1 = sum(value1)/nop;

value2 = mask.*Image2;
value2 = value2(value2~=0);
value2 = reshape(value2,length(value2)/nof,nof);
std2 = std(value2);
value2 = sum(value2)/nop;

fprintf(['mean 1 = ', num2str(mean(value1))])
fprintf(['\nstd 1  = ', num2str(mean(std1))])
fprintf(['\nmean 2 = ', num2str(mean(value2))])
fprintf(['\nstd 2  = ', num2str(mean(std2)),'\n'])
%value2 = value2/mean(value2(:))*mean(value1(:));
figure
imagesc([Image1(:,:,1) + value1(1)*mask,Image2(:,:,1) + value2(1)*mask])
axis image

figure
plot(std1./value1,'LineWidth',2,'Marker','.','MarkerSize',15)
hold on
plot(std2./value2,'LineWidth',2,'Marker','.','MarkerSize',15)
title 'Normalized Standard Deviation'
xlabel 'Time Frame'
set(gca,'FontSize',20)
Axis = axis;
Axis([1,2]) = [1,nof];
axis(Axis);
legend(['Image1 mean=',num2str(mean(std1),'%.2f')],['Image2 mean=',num2str(mean(std2),'%.2f')])
