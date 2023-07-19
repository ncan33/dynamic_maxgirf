function [Cost_new, Cost, fNorm, tNorm] = Cost_TFD_ADMM(fUpdate, Image, l2ref, l2weight, Cost_old)

N = numel(Image);

fNorm = sum(abs(fUpdate(:)).^2);

if l2weight ~= 0
    tNorm = TFD(Image);
    tNorm = tNorm - l2ref;
    tNorm = l2weight * sum(abs(tNorm(:)).^2);
else
    tNorm = 0;
end

fNorm = fNorm/N;
tNorm = tNorm/N;

Cost = tNorm + fNorm;

if nargin == 4
    Cost_new = Cost;
    return
end

Cost_new = Cost_old;

if isempty(Cost_old.fidelityNorm)==1
    Cost_new.fidelityNorm = gather(fNorm);
    Cost_new.temporalNorm = gather(tNorm);
    Cost_new.totalCost = gather(Cost);
else    
    Cost_new.fidelityNorm(end+1) = gather(fNorm);
    Cost_new.temporalNorm(end+1) = gather(tNorm);
    Cost_new.totalCost(end+1) = gather(Cost);
end

end