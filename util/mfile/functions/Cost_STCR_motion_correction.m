function [Cost_new,Cost,fNorm,tNorm,sNorm] = Cost_STCR_motion_correction(fUpdate, Image, Motion, sWeight, tWeight, Cost_old)

%N = numel(Image);

fNorm = sum(abs(fUpdate(:)).^2);

if tWeight ~= 0
    difft = diff(Image,1,3);
    ROI_diff = zeros(size(difft),'like',difft);
    for i = Motion.slice_pick
        temp = Image(:,:,:,:,i);
        temp_shift = Motion.y_motion(:,i);
        for t = 1:size(Image,3)
            temp(:,:,t) = circshift(temp(:,:,t),-temp_shift(t),1);
        end
        temp = abs(diff(temp,1,3));
        for t=1:size(Image,3)-1
            temp(:,:,t) = circshift(temp(:,:,t),temp_shift(t+1),1);
        end
        ROI_diff(:,:,:,:,i) = temp;
    end
    difft(Motion.Mask(:,:,2:end,:,:)) = ROI_diff(Motion.Mask(:,:,2:end,:,:));
    difft = abs(difft);
    tNorm = tWeight * sum(difft(:));
else
    tNorm = 0;
end

if sWeight ~= 0
    sx_norm = abs(diff(Image,1,2));
    sx_norm(:,end+1,:,:,:)=0;
    sy_norm = abs(diff(Image,1,1));
    sy_norm(end+1,:,:,:,:)=0;
    sNorm = sqrt(abs(sx_norm).^2+abs(sy_norm).^2);
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