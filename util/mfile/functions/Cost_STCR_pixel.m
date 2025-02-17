function [Cost_new,Cost,fNorm,tNorm,sNorm] = Cost_STCR_pixel(fUpdate, Image, Motion, sWeight, tWeight, Cost_old)

%N = numel(Image);

fNorm = sum(abs(fUpdate(:)).^2);

if tWeight ~= 0
    tNorm = tWeight .* abs(crop_half_FOV(Image(Motion.idx_b) - Image(:,:,1:end-1,:,:)));
    %tNorm = abs(diff(crop_half_FOV(Image),1,3));
    tNorm = sum(tNorm(:));
else
    tNorm = 0;
end

Image = crop_half_FOV(Image);

if sWeight ~= 0
    sx_norm = abs(diff(Image,1,2));
    sx_norm(:,end+1,:,:,:)=0;
    sy_norm = abs(diff(Image,1,1));
    sy_norm(end+1,:,:,:,:)=0;
    sNorm = sWeight .* sqrt(abs(sx_norm).^2+abs(sy_norm).^2);
    sNorm = sum(sNorm(:));
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