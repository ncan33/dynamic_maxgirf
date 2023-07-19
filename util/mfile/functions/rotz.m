function r = rotz(alpha)

alpha = alpha/180*pi;
r = [cos(alpha),-sin(alpha),0; sin(alpha),cos(alpha),0; 0,0,1];