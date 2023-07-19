
function [Dic] = Bloch_SASHA(SRTs,T1,flip_angle,nor,TR,PD_flip_angle)

M0 = 1;
Mi = 0;
for iSRT=1:length(SRTs)
    Mz = M0 + (Mi-M0)*exp(-(1:SRTs(iSRT))./T1');
    Signal = zeros(length(T1),1);
    for i=1:nor
        Mz(:,end+1) = Mz(:,end).*cos(flip_angle/180*pi);
        Signal(:,end+1) = Mz(:,end-1).*sin(flip_angle/180*pi);
        Mz(:,end+1) = M0 + (Mz(:,end)-M0).*exp(-(TR)./T1');
    end
    Signal(:,1) = [];
    Dic(:,iSRT) = mean(Signal,2);
end


M0 = 1;
Mi = 1;

Mz = M0 + (Mi-M0)*exp(-0./T1');
Signal = zeros(length(T1),1);
for i=1:nor
    Mz(:,end+1) = Mz(:,end).*cos(PD_flip_angle/180*pi);
    Signal(:,end+1) = Mz(:,end-1).*sin(PD_flip_angle/180*pi);
    Mz(:,end+1) = M0 + (Mz(:,end)-M0).*exp(-(TR)./T1');
end
Signal(:,1) = [];
Dic_PD = mean(Signal,2);

Dic = Dic./Dic_PD;


end