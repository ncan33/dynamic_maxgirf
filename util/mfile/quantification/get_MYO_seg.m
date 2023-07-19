function [mask_seg,tissue_curves] = get_MYO_seg(Image,Number_of_seg)

if ~exist('Number_of_seg')
    Number_of_seg = 6;
end
Image = abs(squeeze(Image));

Image_std = mean(Image,3);
fprintf('Draw the outer myocardium \n')
Mask_out = get_mask(Image_std);
fprintf('Draw the inner myocardium \n')
Mask_in = get_mask(Image_std);

Mask_MYO = Mask_out & ~Mask_in;

figure(1)
clf
imagesc(Image_std + Mask_MYO*max(Image_std(:)));
axis image

%% get center of weight of Mask_MYO
center = regionprops(Mask_MYO,'Centroid');
center = center.Centroid;
hold on
plot(center(1),center(2),'r*')
%% get reference point
[x,y] = ginput(1);

r = sqrt(sum((center - [x,y]).^2));
r = r*2;

angle_ref = atan2(center(2)-y,center(1)-x);
angle_all = (0:2*pi/Number_of_seg:2*pi-2*pi/Number_of_seg) + angle_ref;

x_all = center(1) + cos(angle_all)*r;
y_all = center(2) + sin(angle_all)*r;
for i=1:Number_of_seg
    plot([center(1),x_all(i)],[center(2),y_all(i)],'r')
end

angle_all = mod(angle_all,2*pi);
angle_all = angle_all - pi;
[y,x] = find(Mask_MYO);
[theta,r] = cart2pol(x-center(1),y-center(2));

angle_all = sort(angle_all);
mask_seg = false([size(Mask_MYO),Number_of_seg]);
for seg=1:Number_of_seg-1
    mask_temp = theta>angle_all(seg) & theta<angle_all(seg+1);
    theta_temp = theta(mask_temp);
    r_temp = r(mask_temp);
    [x,y] = pol2cart(theta_temp,r_temp);
    y = round(y + center(2));
    x = round(x + center(1));
    mask1 = false(size(Mask_MYO));
    mask1(sub2ind(size(Mask_MYO),y,x)) = true;
    mask_seg(:,:,seg) = mask1;
end
mask_seg(:,:,Number_of_seg) = Mask_MYO & ~(sum(mask_seg,3));

nop = sum(sum(mask_seg));
tissue_curves = sum(sum(mask_seg.*permute(Image,[1,2,4,3])))./nop;
tissue_curves = squeeze(tissue_curves)';
figure
plot(tissue_curves)


d_SI = tissue_curves - mean(tissue_curves(2:10,:));
figure
plot(d_SI)



[~,peak_enhance] = max(sum(tissue_curves,2));
Image_draw = Image(:,:,peak_enhance);

Image_draw = Image_draw./max(Image_draw(:));
Image_draw = Image_draw*1.5;
Image_draw = repmat(Image_draw,[1,1,3]);
%Image_draw = cat(3,Image_draw,zeros(144,144,2));

for i=1:Number_of_seg
    edges(:,:,i) = bwdist(~mask_seg(:,:,i)) == 1;
end

Colors = colors;

for i=1:Number_of_seg
    for j=1:3
        Image_temp = Image_draw(:,:,j);
        Image_temp(edges(:,:,i)) = Colors(i,j);
        Image_draw(:,:,j) = Image_temp;
    end
end
figure
imagesc(Image_draw)
axis image

