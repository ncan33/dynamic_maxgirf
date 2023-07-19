function Image = interp_2_slices_in_between(Image)
% Image should have dimension: x,y,z,frames
[sx,sy,sz,nof] = size(Image);
if nof>1
    [x,y,z,f] = ndgrid(1:sx,1:sy,1:sz,1:nof);
    [xv,yv,zv,fv] = ndgrid(1:sx,1:sy,1:1/3:sz,1:nof);
    Image = interpn(x,y,z,f,Image,xv,yv,zv,fv);
else
    [x,y,z] = ndgrid(1:sx,1:sy,1:sz);
    [xv,yv,zv] = ndgrid(1:sx,1:sy,1:1/3:sz);
    Image = interpn(x,y,z,Image,xv,yv,zv);
end