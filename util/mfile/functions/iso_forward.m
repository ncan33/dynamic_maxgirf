function Image = iso_forward(Image,iso)
[sx,sy,nof,Nslice] = size(Image);
Image = Image(:,:,:,1:3:end);
Image = Image(:,:,:,iso.order_back);
Image = reshape(Image,[sx,sy,nof,1,iso.nSMS,iso.nset]);
return

[sx,sy,nof,Nslice] = size(Image);
Image = reshape(Image,[sx*sy*nof,Nslice]);
% Image = Image*iso.slice_profile;
Nphase = iso.Nphase;
Ncycle = iso.Ncycle;
Image = reshape(Image,[sx,sy,Nphase,Ncycle,Nslice]);
[x,y,ph,cy,z] = ndgrid(1:sx,1:sy,1:Nphase,1:Ncycle,1:Nslice);
[xv,yv,phv,cyv,zv] = ndgrid(1:sx,1:sy,1:Nphase,1:Ncycle,1:(Nslice-1)/(iso.nSMS*iso.nset-1):Nslice);
Image = interpn(x,y,ph,cy,z,Image,xv,yv,phv,cyv,zv);
Image = reshape(Image,[sx,sy,nof,iso.nSMS*iso.nset]);
Image = Image(:,:,:,iso.order_back);
Image = reshape(Image,[sx,sy,nof,1,iso.nSMS,iso.nset]);
