function [Cost_new,Cost,fNorm,tNorm,sNorm] = Cost_STCR_exact(Image, Data, sWeight, tWeight, Cost_old)

N = numel(Image);

fNorm = fft2(Image).*Data.mask - Data.kSpace;
fNorm = sum(abs(fNorm(:)).^2)/N/288^2;

if tWeight ~= 0
    tNorm = abs(diff(Image,1,3));
    tNorm = tWeight * sum(tNorm(:).^2)/N;
else
    tNorm = 0;
end

if sWeight ~= 0
    sx_norm = abs(diff(Image,1,2));
    sx_norm(:,end+1,:,:,:)=0;
    sy_norm = abs(diff(Image,1,1));
    sy_norm(end+1,:,:,:,:)=0;
    sNorm = sqrt(abs(sx_norm).^2+abs(sy_norm).^2);
    sNorm = sWeight * sum(sNorm(:))/N;
else
    sNorm = 0;
end

Cost = sNorm + tNorm + fNorm;

if nargin == 4
    Cost_new = Cost;
    return
end

Cost_new = Cost_old;

if isempty(Cost_old.fidelityNorm)==1
    Cost_new.fidelityNorm = gather(fNorm);
    Cost_new.temporalNorm = gather(tNorm);
    Cost_new.spatialNorm = gather(sNorm);
    Cost_new.totalCost = gather(Cost);
else    
    Cost_new.fidelityNorm(end+1) = gather(fNorm);
    Cost_new.temporalNorm(end+1) = gather(tNorm);
    Cost_new.spatialNorm(end+1) = gather(sNorm);
    Cost_new.totalCost(end+1) = gather(Cost);
end

end