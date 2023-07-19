function dr = Radial_Derivative(Image,center,points)

[sx,sy] = size(Image);

CartCoor = points - center;
[theta,r] = cart2pol(CartCoor(:,1),CartCoor(:,2));
[CartCoor_outerx,CartCoor_outery] = pol2cart(theta,r+sqrt(2));
CartCoor_outerx = round(CartCoor_outerx + center(1));
CartCoor_outery = round(CartCoor_outery + center(2));
CartCoor_outerx(CartCoor_outerx>sx) = sx;
CartCoor_outery(CartCoor_outery>sy) = sy;
CartCoor_outer = [CartCoor_outerx,CartCoor_outery];
CartCoor_outer(CartCoor_outer<=0) = 1;

%idx = abs(points - CartCoor_outer);
%idx = sum(idx,2);
%idx = find(idx==0);
inner_value = Image(sub2ind([sx,sy],points(:,1),points(:,2)));
outer_value = Image(sub2ind([sx,sy],CartCoor_outer(:,1),CartCoor_outer(:,2)));
dr = (inner_value - outer_value)./(inner_value + outer_value);