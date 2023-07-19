
function [signal_1,S1] = Bloch_IR_continupous(InversionTime,InversionEfficiency,T1,siz,flip_angle,nor,TR)

M0 = 1;
Mi = -M0*InversionEfficiency;

% TR = 2.79;
% TE = 1.1210;
% 
% nor = 150;

% if InversionTime(1) < 1
%     InversionTime = InversionTime*1000;
% end

% recovery before measurement 1
S1 = M0 + (Mi-M0)*exp(-(InversionTime(1))/T1);

% measurement 1
for i=1:nor
    S1(end+1) = S1(end)*cos(flip_angle/180*pi);
    %S1(end+1) = M0 + (S1(end)-M0)*exp(-(TE)/T1);
    S1(end+1) = M0 + (S1(end)-M0)*exp(-(TR)/T1);
end
s_temp_1 = S1((end-nor*2):2:end-1);
s_temp_1 = reshape(s_temp_1,siz);
s_temp_1 = abs(mean(s_temp_1,1));
signal_1 = s_temp_1;
