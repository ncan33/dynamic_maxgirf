function [Data,para] = get_iso_image(Data,para)
% slice_thickness = 7;
% inplane_reso = 1.8;

if isfield(Data,'first_guess')
    Image = Data.first_guess;
else
    Image = Data.first_est;
end

[sx,sy,nof,~,nSMS,nset] = size(Image);
nslice = nSMS*nset;

Image = reshape(Image,[sx,sy,nof,nslice]);

order = 1:nSMS:nslice;
for i=nSMS:-1:2
    order = [order,i:nSMS:nslice];
end
[~,order_back] = sort(order);

Image = Image(:,:,:,order);

Nphase = size(para.Recon.bins,1);
Ncycle = nof/Nphase;

Image = reshape(Image,[sx,sy,Nphase,Ncycle,nslice]);

[x,y,ph,cy,z] = ndgrid(1:sx,1:sy,1:Nphase,1:Ncycle,1:nslice);
%[xv,yv,phv,cyv,zv] = ndgrid(1:sx,1:sy,1:Nphase,1:Ncycle,1:(iso.nSMS*iso.nset-1)/(Nslice-1):iso.nSMS*iso.nset);
[xv,yv,phv,cyv,zv] = ndgrid(1:sx,1:sy,1:Nphase,1:Ncycle,1:1/3:nslice);
Image = interpn(x,y,ph,cy,z,Image,xv,yv,phv,cyv,zv);

Nslice = size(Image,5);
%{
load('/v/raid1b/ytian/MRIdata/TestData/190214_Test_isotropic/SliceProfile/SliceNormal.mat')

slice_profile = Slice.Profile;
slice_profile_all = zeros(Nslice,nslice);
slice_profile = [slice_profile;zeros((nslice-5)*200,1)];
N_sp_interp = ceil((nslice*200/Nslice))*Nslice;
interp_coor = 1:(nslice*200-1)/(N_sp_interp-1):nslice*200;
for i=1:nslice
    shift_n = -2+i-1;
    slice_profile_temp = circshift(slice_profile,[shift_n*200,0]);
    if i<3
        slice_profile_temp((5+shift_n)*200+1:end) = 0;
    elseif i>nslice-2
        slice_profile_temp(1:(shift_n+5-nslice)*200) = 0;
    end
    slice_profile_temp = interp1(1:1800,slice_profile_temp,interp_coor);
    slice_profile_temp = sum(reshape(slice_profile_temp,N_sp_interp/Nslice,Nslice),1);
    slice_profile_temp = slice_profile_temp/sum(slice_profile_temp);
    slice_profile_all(:,i) = slice_profile_temp;
end
%}
Image = reshape(Image,[sx,sy,nof,Nslice]);
Data.first_guess = Image;
%Data.iso.slice_profile = slice_profile_all;
%Data.iso.slice_profile_pinv = pinv(slice_profile_all);
Data.iso.Nphase = Nphase;
Data.iso.Ncycle = Ncycle;
Data.iso.order = order;
Data.iso.order_back = order_back;
Data.iso.nSMS = nSMS;
Data.iso.nset = nset;