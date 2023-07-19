function [signal_1,Mz] = Bloch_2_IR_1_SLG(InversionTime,A,T1,siz,flip_angle,SecondInversionTime,nor,TR,TE)

M0 = 1;
Mi = -M0*A;

if InversionTime(1) < 1
    InversionTime = InversionTime*1000;
end

% recovery before measurement 1
Mz = M0 + (Mi-M0)*exp(-(1:InversionTime(1))/T1);

% measurement 1
for i=1:nor
    Mz(end+1) = Mz(end)*cos(flip_angle/180*pi);
    %S1(end+1) = M0 + (S1(end)-M0)*exp(-(TE)/T1);
    Mz(end+1) = M0 + (Mz(end)-M0)*exp(-(TR)/T1);
end
s_temp_1 = Mz((end-nor*2):2:end-1);
s_temp_1 = reshape(s_temp_1,siz);
s_temp_1 = abs(mean(s_temp_1,1));
signal_1 = s_temp_1;

% measurement 2-5
for InvTime = 2:5
    t = 1:(InversionTime(InvTime)-InversionTime(InvTime-1)-nor*TR);
    Mz(end+1:end+length(t)) = M0 + (Mz(end)-M0)*exp(-(t)/T1);
    for i=1:nor
        Mz(end+1) = Mz(end)*cos(flip_angle/180*pi);
        %S1(end+1) = M0 + (S1(end)-M0)*exp(-(TE)/T1);
        Mz(end+1) = M0 + (Mz(end)-M0)*exp(-(TR)/T1);
    end
    s_temp_1 = Mz((end-nor*2):2:end-1);
    s_temp_1 = reshape(s_temp_1,siz);
    s_temp_1 = abs(mean(s_temp_1,1));
    signal_1(end+1:end+siz(2)) = s_temp_1;
end

% Recovery from 5th measuremetn to 2nd inversion
t = 1:SecondInversionTime-InversionTime(5)-nor*TR;
Mz(end+1:end+length(t)) = M0 + (Mz(end)-M0)*exp(-(t)/T1);

% second inversion
Mz(end+1) = - Mz(end)*A;

% recovery between 2nd IR and measurement 6
t = 1:InversionTime(6) - SecondInversionTime;
Mz(end+1:end+length(t)) = M0 + (Mz(end)-M0)*exp(-(t)/T1);

% measurement 6
for i=1:nor
    Mz(end+1) = Mz(end)*cos(flip_angle/180*pi);
    %S1(end+1) = M0 + (S1(end)-M0)*exp(-(TE)/T1);
    Mz(end+1) = M0 + (Mz(end)-M0)*exp(-(TR)/T1);
end
s_temp_1 = Mz((end-nor*2):2:end-1);
s_temp_1 = reshape(s_temp_1,siz);
s_temp_1 = abs(mean(s_temp_1,1));
signal_1(end+1:end+length(s_temp_1)) = s_temp_1;

% measurement 7-8
for InvTime = 7:8
    t = 1:(InversionTime(InvTime)-InversionTime(InvTime-1)-nor*TR);
    Mz(end+1:end+length(t)) = M0 + (Mz(end)-M0)*exp(-(t)/T1);

    for i=1:nor
        Mz(end+1) = Mz(end)*cos(flip_angle/180*pi);
        %S1(end+1) = M0 + (S1(end)-M0)*exp(-(TE)/T1);
        Mz(end+1) = M0 + (Mz(end)-M0)*exp(-(TR)/T1);
    end
    s_temp_1 = Mz((end-nor*2):2:end-1);
    s_temp_1 = reshape(s_temp_1,siz);
    s_temp_1 = abs(mean(s_temp_1,1));
    signal_1(end+1:end+siz(2)) = s_temp_1;
end

signal_1 = signal_1*A;

end