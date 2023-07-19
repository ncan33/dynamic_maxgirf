function [Cost_new,Cost,fNorm,tNorm,sNorm] = Cost_STCR_pixel_dstv(fUpdate, Image, Data, sWeight, tWeight, Cost_old)

%N = numel(Image);

fNorm = sum(abs(fUpdate(:)).^2);

if tWeight ~= 0
    tNorm = abs(crop_half_FOV(Image(Data.Motion.idx_b) - Image(:,:,1:end-1,:,:)));
    %tNorm = abs(diff(crop_half_FOV(Image),1,3));
    tNorm = tWeight * sum(tNorm(:));
else
    tNorm = 0;
end

Image = crop_half_FOV(Image);
edge = crop_half_FOV(Data.edge);
phi = crop_half_FOV(Data.phi);
if sWeight ~= 0
    sx_norm = abs(diff(Image,1,2));
    sx_norm(:,end+1,:,:,:)=0;
    sy_norm = abs(diff(Image,1,1));
    sy_norm(end+1,:,:,:,:)=0;
    sNorm = sqrt(abs(sx_norm).^2+abs(sy_norm).^2).*edge;
    sNorm = sNorm + sqrt(abs(sx_norm.*cos(phi)).^2 + abs(sy_norm.*sin(phi)).^2).*(1-edge);
    sNorm = sWeight * sum(sNorm(:));
else
    sNorm = 0;
end

Cost = sNorm + tNorm + fNorm;

if nargin == 5
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