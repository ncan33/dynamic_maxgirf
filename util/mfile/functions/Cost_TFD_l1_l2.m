function [Cost_new,Cost,fNorm,tNorm,sNorm] = Cost_TFD_l1_l2(fUpdate, Image, l2Weight, tWeight, Cost_old)

N = numel(Image);

fNorm = sum(abs(fUpdate(:)).^2);
% Image = crop_half_FOV(Image);
if l2Weight ~= 0
    l2Norm = mean(l2Weight(:)) .* abs(diff(Image,1,3)).^2;
    l2Norm = sum(l2Norm(:));
else
    l2Norm = 0;
end

if tWeight ~= 0
    tNorm = mean(tWeight(:)) .* abs(diff(Image,1,3));
    tNorm = sum(tNorm(:));
else
    tNorm = 0;
end

% if sWeight ~= 0
%     sx_norm = abs(diff(Image,1,2));
%     sx_norm(:,end+1,:,:,:)=0;
%     sy_norm = abs(diff(Image,1,1));
%     sy_norm(end+1,:,:,:,:)=0;
%     sNorm = sWeight .* sqrt(abs(sx_norm).^2+abs(sy_norm).^2);
%     sNorm = sum(sNorm(:));
% else
%     sNorm = 0;
% end

fNorm = fNorm/N;
tNorm = tNorm/N;
l2Norm = l2Norm/N;

Cost = l2Norm + tNorm + fNorm;

if nargin == 4
    Cost_new = Cost;
    return
end

Cost_new = Cost_old;

if isempty(Cost_old.fidelityNorm)==1
    Cost_new.fidelityNorm = gather(fNorm);
    Cost_new.temporalNorm = gather(tNorm);
    Cost_new.l2Norm = gather(l2Norm);
    Cost_new.totalCost = gather(Cost);
else    
    Cost_new.fidelityNorm(end+1) = gather(fNorm);
    Cost_new.temporalNorm(end+1) = gather(tNorm);
    Cost_new.l2Norm(end+1) = gather(l2Norm);
    Cost_new.totalCost(end+1) = gather(Cost);
end

end