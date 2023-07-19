function PolarImage = GetPolarImage(Image,x_0,y_0,R,R_n,Theta_n)

% Image: input rectangular image
% x_0,y_0: where to start do polar
% R: radius in oringional coordinate
% R_n: number of points along each ray
% Theta_n: number of angle views

[sx,sy] = size(Image);
PolarImage = zeros(R_n,Theta_n);

dR = R/R_n;
dTheta = 2*pi/Theta_n;

for Nr = 1:R_n
for Nt = 1:Theta_n

r = (Nr-1)*dR;
theta = (Nt-1)*dTheta;
dx = r*cos(theta);
dy = r*sin(theta);
x = x_0+dx;
y = y_0+dy;

if y>sy || x>sx || x<1 || y<1
PolarImage(Nr,Nt) = 0;
else
PolarImage(Nr,Nt) = interpolate(Image,x,y);
end

end
end
end

function v = interpolate (imR, x, y)
    xf = floor(x);
    xc = ceil(x);
    yf = floor(y);
    yc = ceil(y);
    if xf == xc && yc == yf
        v = imR (xc, yc);
    elseif xf == xc
        v = imR (xf, yf) + (y - yf)*(imR (xf, yc) - imR (xf, yf));
    elseif yf == yc
        v = imR (xf, yf) + (x - xf)*(imR (xc, yf) - imR (xf, yf));
    else
       A = [ xf yf xf*yf 1
             xf yc xf*yc 1
             xc yf xc*yf 1
             xc yc xc*yc 1 ];
       r = [ imR(xf, yf)
             imR(xf, yc)
             imR(xc, yf)
             imR(xc, yc) ];
       a = A\double(r);
       w = [x y x*y 1];
       v = w*a;
    end
end