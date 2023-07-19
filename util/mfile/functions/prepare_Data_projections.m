function [Data,para] = prepare_Data_projections(kSpace_all,RayPosition,para)

t1 = tic;
[sx,~,no_comp,~,ns] = size(kSpace_all);
%nof                 = para.Recon.nof;
nSMS                = para.Recon.nSMS;
selected_rays_start = para.selected_rays_start;
selected_rays_end   = para.selected_rays_end;
%ifNUFFT             = para.ifNUFFT;
kCenter             = para.kSpace_center;
%interp_method       = para.Recon.interp_method;

nor = selected_rays_end - selected_rays_start + 1; para.Recon.nor = nor;

if para.cSMS
    temp = find(RayPosition==1);
    RayPosition(temp(end-2):end) = 0; % delete the last few rays
    RayPosition(temp(1):temp(2)-1) = 0; %deleta first few rays
    
    index_1 = find(RayPosition==1);
    nor_all = diff(index_1);
    nof_low_ray_number = find(nor_all<selected_rays_end);
    
    for i=1:length(nof_low_ray_number)
        RayPosition(index_1(nof_low_ray_number(i)):index_1(nof_low_ray_number(i)+1)-1) = 0;
    end
end

RayIndex = logical((RayPosition > selected_rays_start-1) .* (RayPosition < selected_rays_end+1));
theta_all     = get_angle_mod(para);% radial sampling angle
phase_mod_all = single(get_phase_mod(para));% if SMS

theta = theta_all(RayIndex); nof = length(theta)/nor;
theta = reshape(theta,[1,nor,nof]);
kSpace_radial = kSpace_all(:,RayIndex,:,:,:);
kSpace_radial = reshape(kSpace_radial,[sx,nor,nof,no_comp,1,ns]);
phase_mod = phase_mod_all(1,RayIndex,1,:);
phase_mod = reshape(phase_mod,[1,nor,nof,1,nSMS]);

if ~isfield(para.Recon,'RF_frames') && isfield(para.Recon,'PD_frames')
    para.Recon.RF_frames = para.Recon.PD_frames(end)+1:nof;
end

%%%%% pre-interpolation
disp('Pre-interpolate onto Cartesian space...')

matObj = matfile(para.dir.PCA_dir);% detect if there's trajectory
varlist = who(matObj,'Kx');
if para.asymmetry
    sx = para.AsymmetryRayEnd;
end
if ~isempty(varlist)
    load(para.dir.PCA_dir,'Kx','Ky')
    kx = double(squeeze(Kx(:,RayIndex))); ky = double(squeeze(Ky(:,RayIndex)));
    kx = reshape(kx,[sx,nor,nof]); ky = reshape(ky,[sx,nor,nof]);
else
    [kx, ky] = get_k_coor(sx,theta,0,kCenter);
end

varlist = who(matObj,'kSpace_info');
if ~isempty(varlist)
    load(para.dir.PCA_dir,'kSpace_info')
    if isfield(kSpace_info,'kx')
        kx = double(squeeze(kSpace_info.kx(:,RayIndex))); ky = double(squeeze(kSpace_info.ky(:,RayIndex)));
        kx = reshape(kx,[sx,nor,nof]); ky = reshape(ky,[sx,nor,nof]);
    end
else
    [kx, ky] = get_k_coor(sx,theta,0,kCenter);
end

[Data.G,Data.kSpace] = GROG_seperate_SMS_GNUFFT(squeeze(kSpace_radial),kx,ky,para);

if nSMS == 1 
    para.Recon.type = '2D';
else
    para.Recon.type = 'seperate SMS';
end

s = round(Data.G.core_size/2);
Data.kSpace = Data.kSpace(s+1:s+Data.G.sx_over, s+1:s+Data.G.sx_over,:,:,:,:,:);

if para.image_orintation == 0
    para.image_orintation = orintation_detection(abs(fftshift(ifft2(sum(sum(sum(Data.kSpace,3),4),7)))));
    Data.kSpace = orintate_image(Data.kSpace,para.image_orintation);
else
    Data.kSpace = orintate_image(Data.kSpace,para.image_orintation);
end

Data.kSpace = fftshift2(Data.kSpace);
Data.mask = logical(abs(Data.kSpace(:,:,:,1,:,:,:)));

im = ifft2(Data.kSpace);
im = fftshift2(im);
Data.kSpace = fft2(im).*Data.mask;
sx_over = size(im,1);
para.Recon.kSpace_size = [sx_over,sx_over];
Data.sens_map = get_sens_map(im,'2D');

para.Recon.nof = para.Recon.nof*para.Recon.nor;
para.Recon.nor = 1;

Data.filter = ramp_filter_for_pre_interp(para);
        
im = ifft2(Data.kSpace.*Data.filter);

Data.first_est = single(sum(bsxfun(@times, im, conj(Data.sens_map)),4));

scale_image = mean(abs(Data.first_est(:)));
para.Recon.weight_tTV = scale_image*para.weight_tTV;
para.Recon.weight_sTV = scale_image*para.weight_sTV;

para.CPUtime.prepare_kSpace = toc(t1);toc(t1);fprintf('\n');

function [G,kSpace_cart] = GROG_seperate_SMS_GNUFFT(kSpace_radial, kx, ky, para)

[sx,nor,nof,nc,~,NSlice] = size(kSpace_radial);
core_size = para.core_size;
over_sampling = para.over_sampling;

for i=1:NSlice
    kSpace_radial_temp = kSpace_radial(:,:,:,:,i);
    kSpace_radial_temp = reshape(kSpace_radial_temp,[sx nor nof nc]);
    kx = reshape(kx,[sx nor nof]);
    ky = reshape(ky,[sx nor nof]);
    G = GNUFFT_init(kSpace_radial_temp,kx,ky,over_sampling,core_size);
    kSpace_radial_temp = reshape(kSpace_radial_temp,[sx,1,nor*nof,nc]);
    kSpace_cart(:,:,:,:,:,i) = GROG.GNUFFT_rad2cart(kSpace_radial_temp,G);
end

sx_over = size(kSpace_cart,1);
kSpace_cart = reshape(kSpace_cart,[sx_over,sx_over,nof*nor,nc,1,NSlice]);

function G = GNUFFT_init(kSpace_radial,kx,ky,over_sampling,core_size)

sx  = size(kSpace_radial,1);
nx = round(sx/8);
[Gx,Gy] = GROG.get_Gx_Gy(kSpace_radial(nx*3+1:5*nx,:,:,:), kx(nx*3+1:5*nx,:,:,:), ky(nx*3+1:5*nx,:,:,:));

skx = size(kx,1);
nor = size(kx,2);
nof = size(kx,3);
nc  = size(Gx,1);

%kSpace_radial = reshape(kSpace_radial,[skx,1,nor*nof,nc]);
kx = reshape(kx,[skx,1,nor*nof]);
ky = reshape(ky,[skx,1,nor*nof]);
nof = nor*nof;
nor = 1;

sx_over = round(sx*over_sampling);

kx = kx*over_sampling;
ky = ky*over_sampling;
Gx = Gx^(1/over_sampling);
Gy = Gy^(1/over_sampling);

max_shift_x = core_size(1)/2;
max_shift_y = core_size(2)/2;

GxDict_size = max_shift_x*200 + 1;
GyDict_size = max_shift_y*200 + 1;
GxDict = single(zeros([1 nc nc GxDict_size]));
GyDict = single(zeros([nc nc 1 1 GyDict_size]));

dx = -max_shift_x:0.01:max_shift_x;
dy = -max_shift_y:0.01:max_shift_y;

for di=1:GxDict_size
    GxDict(1,:,:,di) = Gx^dx(di);
end
for di=1:GyDict_size
    GyDict(:,:,1,1,di) = Gy^dy(di);
end

G_r2c_Dict = bsxfun(@times,GxDict,GyDict);
G_r2c_Dict = squeeze(sum(G_r2c_Dict,2));
G_r2c_Dict = reshape(G_r2c_Dict,[nc nc GxDict_size*GyDict_size]);
G_r2c_Dict = permute(G_r2c_Dict,[3 1 2]);

GxDict = permute(GxDict,[2 3 1 4]);
GyDict = permute(GyDict,[3 1 2 4 5]);
G_c2r_Dict = bsxfun(@times,GxDict,GyDict);
G_c2r_Dict = squeeze(sum(G_c2r_Dict,2));
G_c2r_Dict = reshape(G_c2r_Dict,[nc nc GxDict_size*GyDict_size]);
G_c2r_Dict = permute(G_c2r_Dict,[3 1 2]);

x_cart = zeros(skx,nor,nof,core_size(1));
y_cart = zeros(skx,nor,nof,core_size(2));

if mod(core_size(1),2) == 0
    x_cart(:,:,:,1) = floor(kx) - core_size(1)/2 + 1;
else
    x_cart(:,:,:,1) = round(kx) - (core_size(1)-1)/2;
end
add_x(1,1,1,:) = 1:core_size(1)-1;
x_cart(:,:,:,2:end) = bsxfun(@plus,x_cart(:,:,:,1),add_x);

if mod(core_size(2),2) == 0
    y_cart(:,:,:,1) = floor(ky) - core_size(2)/2 + 1;
else
    y_cart(:,:,:,1) = round(ky) - (core_size(2)-1)/2;
end
add_y(1,1,1,:) = 1:core_size(2)-1;
y_cart(:,:,:,2:end) = bsxfun(@plus,y_cart(:,:,:,1),add_y);

x_cart = repmat(x_cart,1,1,1,core_size(2));
x_cart = reshape(x_cart,skx,nor,nof,prod(core_size));
y_cart = permute(y_cart,[1,2,3,5,4]);
y_cart = repmat(y_cart,1,1,1,core_size(1));
y_cart = reshape(y_cart,skx,nor,nof,prod(core_size));

dx = round((x_cart - kx)*100)/100;
dy = round((y_cart - ky)*100)/100;

xDict_r2c = round((dx+max_shift_x)*100)+1;
yDict_r2c = round((dy+max_shift_y)*100)+1;

index_Dict_r2c = sub2ind([GxDict_size GyDict_size],xDict_r2c,yDict_r2c);

xDict_c2r = round(-dx+max_shift_x)*100+1;
yDict_c2r = round(-dy+max_shift_y)*100+1;

index_Dict_c2r = sub2ind([GxDict_size GyDict_size],xDict_c2r,yDict_c2r);

if mod(core_size(1),2) == 0
    x_cart = x_cart+sx_over/2+core_size(1)/2;
else
    x_cart = x_cart+sx_over/2+core_size(1)/2+0.5;
end
if mod(core_size(2),2) == 0
    y_cart = y_cart+sx_over/2+core_size(1)/2;
else
    y_cart = y_cart+sx_over/2+core_size(1)/2+0.5;
end

x_cart = mod(x_cart,sx_over); x_cart(x_cart==0) = sx_over;
y_cart = mod(y_cart,sx_over); y_cart(y_cart==0) = sx_over;

if mod(core_size(1),2) == 0
    indx = sub2ind([(sx_over+core_size(1)),(sx_over+core_size(2)),nof],x_cart,y_cart);
else
    indx = sub2ind([(sx_over+core_size(1)+1),(sx_over+core_size(2)+1),nof],x_cart,y_cart);
end

rad2cart = cell(1,nof);
rad_num = (1:skx*nor*prod(core_size)).';

weight_scale = sqrt(max_shift_x^2+max_shift_y^2)*1.00001;
weight = weight_scale - sqrt(dx.^2 + dy.^2);


weight = permute(weight,[1,2,4,3]);
weight = reshape(weight,sx*nor*prod(core_size),nof);


if mod(core_size(1),2) == 0
    weight_cart = single(zeros((sx_over+core_size(1))*(sx_over+core_size(2)),nof));
else
    weight_cart = single(zeros((sx_over+core_size(1)+1)*(sx_over+core_size(2)+1),nof));
end

for i=1:nof
    indx_temp = indx(:,:,i,:);
    indx_temp = indx_temp(:);
    if mod(core_size(1),2) == 0
        rad2cart{i} = sparse(indx_temp,rad_num,1,(sx_over+core_size(1))*(sx_over+core_size(2)),skx*nor*prod(core_size));
    else
        rad2cart{i} = sparse(indx_temp,rad_num,1,(sx_over+core_size(1)+1)*(sx_over+core_size(2)+1),skx*nor*prod(core_size));
    end
    weight_cart(:,i) = rad2cart{i} * double(weight(:,i));
    
    weight_temp = full(sum(rad2cart{i},2));
    weight_temp(weight_temp ==0) = 1;
    weight3(:,i) = 1./weight_temp;
    clear indx_temp
end


G.rad2cart = rad2cart;
G.Dict_r2c = G_r2c_Dict;
G.Dict_c2r = G_c2r_Dict;
G.indx_r2c = index_Dict_r2c;
G.indx_c2r = index_Dict_c2r;
G.core_size = core_size;
G.sx_over = sx_over;
G.siz = [skx nor nof nc];
G.weight = weight3;
G.weight1 = reshape(weight,[skx*nor,prod(core_size),nof]);
G.weight1 = permute(G.weight1,[1 3 2]);
weight_cart(weight_cart==0)=1;
G.weight2 = 1./weight_cart;

%G.W = density_comp_area(kx,ky,mod(atan(ky(1,:,:)./kx(1,:,:)),pi));
%G.W = zeros(sx_over,sx_over,'single');
%G.W(round(sx_over/2),round(sx_over/2)) = 1;
%G.W = fftshift2(bwdist(G.W));
