function [PolarImage,X,Y] = GetPolarImage_new(Image,x_0,y_0,R,R_n,Theta_n)

% Image: input rectangular image
% x_0,y_0: where to start do polar
% R: radius in oringional coordinate
% R_n: number of points along each ray
% Theta_n: number of angle views

[sx,sy] = size(Image);

dR = R/R_n;
dTheta = 2*pi/Theta_n;

[R,Theta] = meshgrid(0:dR:R-dR,0:dTheta:2*pi-dTheta);

X = R.*cos(Theta);
Y = R.*sin(Theta);

X = X + x_0;
Y = Y + y_0;

%OutPoints = X>sx | Y>sy | X<1 | Y<1;
%X(OutPoints) = nan;
%Y(OutPoints) = nan;

[x,y] = meshgrid(1:sx,1:sy);
PolarImage = griddata(x,y,double(Image),X,Y,'nearest');