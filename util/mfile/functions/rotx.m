function r = rotx(alpha)

alpha = alpha/180*pi;
r = [1,0,0; 0,cos(alpha),-sin(alpha); 0,sin(alpha),cos(alpha)];