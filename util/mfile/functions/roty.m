function r = roty(alpha)

alpha = alpha/180*pi;
r = [cos(alpha),0,sin(alpha); 0,1,0; -sin(alpha),0,cos(alpha)];