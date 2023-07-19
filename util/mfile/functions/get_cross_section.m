function out = get_cross_section(Image)

siz = size(Image);
figure
imagesc(sum(Image(:,:,1,:),4))
axis image
colormap gray
brighten(0.4)

h = imline;

coor = getPosition(h);
r = sqrt((coor(1)-coor(2))^2+(coor(3)-coor(4))^2);

x = coor(1):(coor(2)-coor(1))/r:coor(2);
y = coor(3):(coor(4)-coor(3))/r:coor(4);



x = repmat(x.',[1,siz(3)]);
y = repmat(y.',[1,siz(3)]);
z = repmat(1:siz(3),[length(x),1]);




[X,Y,Z] = meshgrid(1:siz(1),1:siz(2),1:siz(3));
for i=1:siz(4)
    out(:,:,i) = interp3(X,Y,Z,Image(:,:,:,i),x,y,z);
end