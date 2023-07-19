function [center,Mask_LV] = get_seed(Image,LV_seed)
n_theta = 180;
n_r = 30;

[sx,sy] = size(Image);
LV = LV_seed;% this value should be from FindLVRV

%[PolarIm,X,Y] = GetPolarImage_new(Image,LV(1),LV(2),n_r,n_r,n_theta);

%threshood = Image(LV(1),LV(2))*0.6;

dist = false(sx,sy);
dist(LV(1),LV(2)) = true;
dist = bwdist(dist);
Mask_LV = seg_LV_yt(Image,dist<15,100,0.5,0,0);
Mask_LV = refine_mask_edge_new(Mask_LV);
%Mask_MY = bwdist(Mask_LV)<6;
%Mask_MY_Rect = seg_MY_yt(Image,Mask_MY,200,10,0.4,0);

s = regionprops(Mask_LV);
center = round(s(1).Centroid);
%{
[PolarIm,X,Y] = GetPolarImage_new(Image,center(1),center(2),n_r,n_r,n_theta);
Mask_Polar_LV = GetPolarImage_new(Mask_LV,center(1),center(2),n_r,n_r,n_theta);
kernel = [1 1 -1 -1];
for i=1:n_theta
    Conv(i,:) = conv(PolarIm(i,:),kernel);
end
Conv(:,1) = [];
Conv(:,end-1:end) = [];
Conv = Conv.*~Mask_Polar_LV;
%x = Boundry_LV(:,2) - center(1);
%y = Boundry_LV(:,1) - center(2);
%r = sqrt(x.^2+y.^2);
%r = median(r);

[PolarIm,X,Y] = GetPolarImage_new(Image,center(1),center(2),n_r,n_r,n_theta);
Mask = PolarIm>threshood;
[Mask,N] = bwlabel(Mask);
for i=1:N
    Area(i) = sum(sum(Mask==i));
end
[~,I] = sort(Area,'descend');
Mask_LV = Mask==Mask(1,1);
Mask_RV = Mask==I(2);
Mask_MY = false(size(Mask));

LV_r_all = sum(Mask_LV,2);
r_polar = median(LV_r_all);
LV_r_all(LV_r_all<r_polar-1 | LV_r_all>r_polar+2) = r_polar;

Mask_LV = false(n_theta,30);

%MY_w = 0;
%MY_c = 0;
for i=1:n_theta
    Mask(i,1:LV_r_all(i)) = 1;
    Mask_LV(i,1:LV_r_all(i)) = true;
    lv = find(Mask_LV(i,:));
    rv = find(Mask_RV(i,:)==true);
    if ~isempty(rv) && ~isempty(lv)
        Mask_MY(i,lv(end)+1:rv(1)-1) = true;
        %MY_w = MY_w + rv(1)-lv(end)-1;
        %MY_c = MY_c+1;
    end
end

%if MY_c~=0
%    MY_w = round(MY_w/MY_c);
%end

MY_r_all = sum(Mask_MY,2);
MY_r_all_nozero = MY_r_all;
MY_r_all_nozero(MY_r_all_nozero==0) = [];
MY_w = median(MY_r_all_nozero);
MY_r_all(MY_r_all<MY_w-1 | MY_r_all>MY_w+1) = MY_w;

Mask_MY = false(n_theta,30);

for i=1:n_theta
    Mask_MY(i,LV_r_all(i)+1:LV_r_all(i)+MY_r_all(i)) = true;
    %if sum(Mask_MY(i,:))==0
    %    temp = Mask(i,:);
    %    lv = find(temp==1);
    %    if lv(end)<r
    %        Mask_MY(i,round(r)+1:round(r)+MY_w) = true;
    %    else
    %        Mask_MY(i,lv(end)+1:lv(end)+MY_w) = true;
    %    end
    %end
end
keyboard
Mask_MY_Rect = seg_LV_yt(Image,Polar2Rect(Mask_LV+Mask_MY,144,144,X,Y),200,10,0.4,0);
Mask_LV_Rect = seg_LV_yt(Image,Polar2Rect(Mask_LV,144,144,X,Y),200,0.2,0.01,0);
Mask_LV = GetPolarImage_new(Mask_LV_Rect,center(1),center(2),n_r,n_r,n_theta);
Mask_MY = GetPolarImage_new(Mask_MY_Rect,center(1),center(2),n_r,n_r,n_theta);

%% get seed

seed_LV_angle = 0:20:359;
N = length(seed_LV_angle);
r = zeros(1,N);
for i=1:N
    for j=1:5
        angle_this = round((seed_LV_angle(i) + j)/2);
        temp = find(Mask_LV(angle_this,:));
        r(i) = r(i)+temp(end);
    end
end
r = r/5;
[seed_lv(1,:),seed_lv(2,:)] = pol2cart(seed_LV_angle/180*pi,r);
seed_lv(1,:) = seed_lv(1,:) + center(1);
seed_lv(2,:) = seed_lv(2,:) + center(2);

seed_MY_angle = seed_LV_angle + round(5);
N = length(seed_MY_angle);
r = zeros(1,N);
for i=1:N
    for j=1:5
        angle_this = round((seed_MY_angle(i) + j)/2);
        temp = find(Mask_MY(angle_this,:));
        r(i) = r(i)+temp(end);
    end
end
r = r/5;
[seed_my(1,:),seed_my(2,:)] = pol2cart(seed_MY_angle/180*pi,r);
seed_my(1,:) = seed_my(1,:) + center(1);
seed_my(2,:) = seed_my(2,:) + center(2);
seed_all = [center.',seed_my,seed_lv];
%}
end