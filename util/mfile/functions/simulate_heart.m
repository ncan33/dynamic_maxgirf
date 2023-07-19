function heart = simulate_heart(siz,cardiac_signal,respiration_signal)

x = siz;
y = siz;

if nargin == 1
    lv_radius = (x+y)/30;
    nof = 1;
else
    nof = length(cardiac_signal);
    lv_radius_dia = (x+y)/15;
    lv_radius_sys = (x+y)/30;
    lv_radius = (lv_radius_dia - lv_radius_sys).*cardiac_signal + lv_radius_sys;
end

lv_blood = false(x,y,nof);
lv_myro  = false(x,y,nof);
rv_blood = false(x,y,nof);

[X,Y] = meshgrid(1:x,1:y);

%im = zeros(x,y);
lv_center_x = round(x/2);
lv_center_y = round(y/2) + round(y/10);
lv_mask = false(x,y);
lv_mask(lv_center_x,lv_center_y) = true;

lv_myo_radius = (x+y)/40+(x+y)/20;
lv_dist = bwdist(lv_mask);


for i=1:nof
    lv_blood(:,:,i) = lv_dist<lv_radius(i);
    lv_myro(:,:,i) = lv_dist>lv_radius(i) & lv_dist<lv_myo_radius;
end

rv_center_x = round(x/2);
rv_center_y = round(y/2) - round(y/15);
%rv_mask = false(x,y);
%rv_mask(rv_center_x,rv_center_y) = true;
rv_a = (x+y)/20;
rv_b = (x+y)/30;
rv_dist = (X-rv_center_y).^2/rv_a.^2 + (Y-rv_center_x).^2/rv_b.^2;
for i=1:nof
    rv_blood(:,:,i) = rv_dist<4;
    rv_blood(:,:,i) = rv_blood(:,:,i) & ~lv_blood(:,:,i);
    rv_blood(:,:,i) = rv_blood(:,:,i) & ~lv_myro(:,:,i);
end


if exist('respiration_signal','var')
    respiration_signal = round((respiration_signal-0.5)*6);
    for i=1:nof
        lv_blood(:,:,i) = circshift(lv_blood(:,:,i),respiration_signal(:,i));
        rv_blood(:,:,i) = circshift(rv_blood(:,:,i),respiration_signal(:,i));
        lv_myro(:,:,i) = circshift(lv_myro(:,:,i),respiration_signal(:,i));
    end
end

heart.lv_blood = lv_blood;
heart.rv_blood = rv_blood;
heart.lv_myro  = lv_myro;
heart.im_size = siz;
heart.lv = [lv_center_x,lv_center_y,lv_myo_radius];
heart.rv = [rv_center_x,rv_center_y,rv_a,rv_b];
heart.lv_radius = lv_radius;

%figure,imagesc(heart.lv_blood*3+heart.lv_myro*2+heart.rv_blood*3);
%colormap gray
%axis image
%axis off