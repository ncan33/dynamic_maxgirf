heart = simulate_heart(128);
[Y,X] = meshgrid(1:heart.im_size,1:heart.im_size);

para.M0     = 1000;
para.alpha  = 12;% flip angle
para.TE     = 1.6;
para.T2star = 50;
para.TR     = 2.7;
para.TD     = 60;

T1_tissue = 1100;
T1_blood  = 1500;
gd_consentration = 4;
gd_relax_tissue = 4;
gd_relax_blood = 5.5;

lv_blood_60_alpha_pulses = Y < heart.lv(2)+heart.lv(3)/2 & Y > heart.lv(2)-heart.lv(3)/2 & heart.lv_blood;
lv_myro_60_alpha_pulses  = Y < heart.lv(2)+heart.lv(3)/2 & Y > heart.lv(2)-heart.lv(3)/2 & heart.lv_myro;

lv_blood_30_alpha_pulses = ~lv_blood_60_alpha_pulses & heart.lv_blood;
lv_myro_30_alpha_pulses = ~lv_myro_60_alpha_pulses & heart.lv_myro;


% pre_contrast
para.T1 = T1_tissue;
para.n  = 15;
lv_myro_30_alpha_pulses = simulate_signal(para)*lv_myro_30_alpha_pulses;

para.n = 30;
lv_myro_60_alpha_pulses = simulate_signal(para)*lv_myro_60_alpha_pulses;

para.T1 = T1_blood;
para.n  = 15;
lv_blood_30_alpha_pulses = simulate_signal(para)*lv_blood_30_alpha_pulses;
rv_blood_30_alpha_pulses = simulate_signal(para)*heart.rv_blood;

para.n = 30;
lv_blood_60_alpha_pulses = simulate_signal(para)*lv_blood_60_alpha_pulses;

im_pre = lv_blood_30_alpha_pulses + lv_blood_60_alpha_pulses + lv_myro_30_alpha_pulses + lv_myro_60_alpha_pulses + rv_blood_30_alpha_pulses;
figure,imagesc(im_pre);

colormap gray
axis image
axis off


% during_contrast

lv_blood_60_alpha_pulses = Y < heart.lv(2)+heart.lv(3)/2 & Y > heart.lv(2)-heart.lv(3)/2 & heart.lv_blood;
lv_myro_60_alpha_pulses  = Y < heart.lv(2)+heart.lv(3)/2 & Y > heart.lv(2)-heart.lv(3)/2 & heart.lv_myro;

lv_blood_30_alpha_pulses = ~lv_blood_60_alpha_pulses & heart.lv_blood;
lv_myro_30_alpha_pulses = ~lv_myro_60_alpha_pulses & heart.lv_myro;

para.T1 = 1/(1/(T1_tissue/1000)+gd_relax_tissue*gd_consentration)*1000;
para.n  = 15;
lv_myro_30_alpha_pulses = simulate_signal(para)*lv_myro_30_alpha_pulses;

para.n = 30;
lv_myro_60_alpha_pulses = simulate_signal(para)*lv_myro_60_alpha_pulses;

para.T1 = 1/(1/(T1_blood/1000)+gd_relax_blood*gd_consentration)*1000;
para.n  = 15;
lv_blood_30_alpha_pulses = simulate_signal(para)*lv_blood_30_alpha_pulses;
rv_blood_30_alpha_pulses = simulate_signal(para)*heart.rv_blood;

para.n = 30;
lv_blood_60_alpha_pulses = simulate_signal(para)*lv_blood_60_alpha_pulses;

im_dur = lv_blood_30_alpha_pulses + lv_blood_60_alpha_pulses + lv_myro_30_alpha_pulses + lv_myro_60_alpha_pulses + rv_blood_30_alpha_pulses;
figure,imagesc(im_dur);

colormap gray
axis image
axis off

% add noise
noise = wgn(heart.im_size,heart.im_size,10);
figure,imagesc(im_pre+noise);

colormap gray
axis image
axis off


figure
imagesc([im_pre,im_dur;im_pre+noise,im_dur+noise])