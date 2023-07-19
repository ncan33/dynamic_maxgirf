function out = simulate_signal(para)

M0     = para.M0;
alpha  = para.alpha/180*pi;% flip angle
TE     = para.TE;
T2star = para.T2star;
T1     = para.T1;
TR     = para.TR;
TD     = para.TD;
n      = para.n; % number of alpha pulses

M = M0*sin(alpha)*exp(-TE/T2star);
a = cos(alpha)*exp(-TR/T1);

out = M*(1-exp(-TD/T1))*a^(n-1)+M*(1-exp(-TR/T1))*(1-a^(n-1))/(1-a);
