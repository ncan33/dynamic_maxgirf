function [mask4] = segment_short_axis_apical(Image)
colors = [0,0.447,0.7410; 0.85,0.325,0.098; 0.929,0.694,0.125; 0.494,0.184,0.556; 0.466,0.674,0.188; 0.301,0.745,0.933];

fprintf('draw LV \n\n')
temp = figure;
set(temp,'Position',[10,10,1000,1000])
imagesc(Image),axis image
h = imfreehand(gca,'closed',false);
lv_mask = createMask(h);

lv_dist = bwdist(lv_mask) - bwdist(~lv_mask);
outside = 6;
inside = 1;
my_mask = bwdist(lv_mask)<outside & bwdist(lv_mask)>inside;


imagesc(Image.*~lv_mask + my_mask*max(max(lv_mask.*Image)))
axis image
colormap jet

clc
prompt = 'Rough adjust contour? \n 1 : more inside \n 2 : less inside \n 3 : more outside \n 4 : less outside \n f : finished \n';
flag = input(prompt,'s');
while ~contains(flag,'f')
    switch flag
        case '1'
            inside = inside - 1;
        case '2'
            inside = inside + 1;
        case '3'
            outside = outside + 1;
        case '4'
            outside = outside - 1;
    end
    my_mask = bwdist(lv_mask)<outside & bwdist(lv_mask)>inside;
    imagesc(Image.*~lv_mask + my_mask*max(max(lv_mask.*Image)))
    axis image
    colormap jet
    flag = input(prompt,'s');
end


clc
prompt = 'Do I need to adjust contour? \n ap : add pixels \n dp : delete pixels \n ar : add area \n dr : delete area \n f : finished \n';
flag = input(prompt,'s');

while ~contains(flag,'f')
    my_mask = adjust_mask(Image,lv_mask,my_mask,flag);
    clc
    flag = input(prompt,'s');
end
b = bwboundaries(my_mask);
imagesc(Image)
axis image
colormap gray
brighten(0.4)
hold on
plot(b{1}(:,2),b{1}(:,1),'LineWidth',2);

center_lv = regionprops(lv_mask,'centroid');
center_lv = center_lv.Centroid;
plot(center_lv(1),center_lv(2),'.','MarkerSize',30)

% for i=1:2
    fprintf('get a reference point (center RV)\n')
    p = impoint;
    p = getPosition(p);
    [theta,r] = cart2pol(p(1)-center_lv(1),p(2)-center_lv(2));
    theta(1) = theta - pi/4;
    theta(2) = theta(1) + pi/2;
    r(2) = r+5;
% end

theta_all = [theta,theta + pi];
r = max(r);

%theta_all = theta + pi/3 *(0:5);
theta_all = mod(theta_all,2*pi) - pi;
[x,y] = pol2cart(theta_all,r+5);

color3 = [0.929,0.694,0.125];
for seg=1:4
    plot([x(seg)+center_lv(1),center_lv(1)],[y(seg)+center_lv(2),center_lv(2)],'color',color3,'LineWidth',2)
end

[y,x] = find(my_mask);
[theta,r] = cart2pol(x-center_lv(1),y-center_lv(2));
%theta = mod(theta);
theta_all = sort(theta_all);
mask4 = false([size(Image),4]);
for seg=1:3
    mask_temp = theta>theta_all(seg) & theta<theta_all(seg+1);
    theta_temp = theta(mask_temp);
    r_temp = r(mask_temp);
    [x,y] = pol2cart(theta_temp,r_temp);
    y = round(y + center_lv(2));
    x = round(x + center_lv(1));
    mask1 = false(size(Image));
    mask1(sub2ind(size(Image),y,x)) = true;
    mask4(:,:,seg) = mask1;
end
mask4(:,:,4) = my_mask & ~(sum(mask4,3));


for seg = 1:4
    [y,x] = find(mask4(:,:,seg));
    plot(x,y,'.','Color',colors(seg,:),'MarkerSize',20)
end
%legend({'1';'2';'3';'4';'5';'6'})

%% LV
% fprintf('Draw LV to caculate AIF: \n')
% h = drawfreehand;
% mask_LV = createMask(h);
%%

figure1(Image,mask4,colors);

clc
flag = input('Do I need to adjust mask4? \n y/n \n','s');
while contains(flag,'y')
    mask4 = adjust_mask4(mask4,Image,colors);
    flag = input('Do I need to adjust mask4? \n y/n \n','s');
end

% tissue_curve = squeeze(mean(mean(current_slice_reg.*permute(mask6,[1,2,4,3]))));
% PD_slice = PD_frame(crop_p(2):crop_p(4),crop_p(1):crop_p(3),:,islice);
% for iframe = 1:size(PD_slice,3)
%     PD_slice(:,:,iframe) = Reg_GS_tv(PD_slice(:,:,iframe),Image,1,50);
% end
% PD_slice = mean(PD_slice,3);
% PD_correction = mean(mean(PD_slice.*mask6));
% PD_correction = permute(PD_correction,[1,3,2]);
% figure,plot(tissue_curve./PD_correction)


function my_mask = adjust_mask(Image,lv_mask,my_mask,flag)

imagesc(Image.*~lv_mask+ my_mask*max(max(lv_mask.*Image)))
colormap jet
axis image

switch flag
    case 'ap'
        [y,x,~] = impixel;
        idx = sub2ind(size(lv_mask),x,y);
        my_mask(idx) = true;
    case 'dp'
        [y,x,~] = impixel;
        idx = sub2ind(size(lv_mask),x,y);
        my_mask(idx) = false;
    case 'ar'
        h = imfreehand;
        mask = createMask(h);
        my_mask = my_mask | mask;
    case 'dr'
        h = imfreehand;
        mask = createMask(h);
        my_mask = my_mask & ~mask;
end
imagesc(Image.*~lv_mask+ my_mask*max(max(lv_mask.*Image)))
colormap jet
axis image

end

function mask4 = adjust_mask4(mask4,Image,colors)

area = input('Which area? \n');
mask_temp = mask4(:,:,area);
prompt = 'How? \n ap : add pixels \n dp : delete pixels \n ar : add area \n dr : delete area \n f : finished \n';
flag = 'y';
clc
%stop = input('Finish changing this area? \n y/n \n','s');
while ~contains(flag,'f')
    flag = input(prompt,'s');
    figure(1)
    switch flag
        case 'ap'
            [y,x,~] = impixel;
            idx = sub2ind(size(mask_temp),x,y);
            mask_temp(idx) = true;
            mask4(:,:,area) = mask_temp;
        case 'dp'
            [y,x,~] = impixel;
            idx = sub2ind(size(mask_temp),x,y);
            mask_temp(idx) = false;
            mask4(:,:,area) = mask_temp;
        case 'ar'
            h = imfreehand;
            mask = createMask(h);
            mask4(:,:,area) = mask4(:,:,area) | mask;
        case 'dr'
            h = imfreehand;
            mask = createMask(h);
            mask4(:,:,area) = mask4(:,:,area) & ~mask;
    end
    mask_temp = mask4(:,:,area);
    figure1(Image,mask4,colors);
    clc
    %flag = input('Finish changing this area? \n y/n \n','s');
end

end


function figure1(Image,mask6,colors)
figure(1)
clf
imagesc(Image),axis image, colormap gray, brighten(0.4),
hold on
for seg = 1:4
    [y,x] = find(mask6(:,:,seg));
    plot(x,y,'.','Color',colors(seg,:),'MarkerSize',15)
end
legend({'1';'2';'3';'4'})
end

end