function Image = iso_backward(Image,iso)

[sx,sy,nof,~,~,~] = size(Image);
%Nslice = size(iso.slice_profile_pinv,2);
Image = reshape(Image,[sx*sy*nof,iso.nSMS*iso.nset]);
Image = Image(:,iso.order);
%Image = Image*iso.slice_profile_pinv;
%Image = reshape(Image,[sx,sy,nof,Nslice]);
Nphase = iso.Nphase;
Ncycle = iso.Ncycle;
Image = reshape(Image,[sx,sy,Nphase,Ncycle,iso.nSMS*iso.nset]);
[x,y,ph,cy,z] = ndgrid(1:sx,1:sy,1:Nphase,1:Ncycle,1:iso.nSMS*iso.nset);
%[xv,yv,phv,cyv,zv] = ndgrid(1:sx,1:sy,1:Nphase,1:Ncycle,1:(iso.nSMS*iso.nset-1)/(Nslice-1):iso.nSMS*iso.nset);
[xv,yv,phv,cyv,zv] = ndgrid(1:sx,1:sy,1:Nphase,1:Ncycle,1:1/3:iso.nSMS*iso.nset);
Image = interpn(x,y,ph,cy,z,Image,xv,yv,phv,cyv,zv);
Image = reshape(Image,[sx,sy,nof,size(Image,5)]);