function [inner_value,outer_value] = get_line_average(Image,idx,center)

[sx,sy] = size(Image);
np = length(idx);
[x,y] = ind2sub([sx,sy],idx);
x = x - center(2);
y = y - center(1);
[theta,r] = cart2pol(x,y);
r_max = round(max(r))+5;
PolarIm = GetPolarImage_new(Image,center(1),center(2),r_max,r_max,180);
theta = round(mod(-theta+pi/2,2*pi)/pi*180/2);
theta(theta==0) = 180;
r = round(r);
r(r<=4) = 5;

for i=1:np
    inner_value(i,1) = mean(PolarIm(theta(i),r(i)-4:r(i)-1));
    outer_value(i,1) = mean(PolarIm(theta(i),r(i)+1:r(i)+4));
end
