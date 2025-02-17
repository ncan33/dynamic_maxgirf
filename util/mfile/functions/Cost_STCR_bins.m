function [Cost_new,Cost,fNorm,tNorm,sNorm] = Cost_STCR_bins(fUpdate, Image, para, Cost_old)

N = numel(Image);
sWeight = para.Recon.weight_sTV;
tWeight = para.Recon.weight_tTV;
bins = para.Recon.bins;

fNorm = sum(abs(fUpdate(:)).^2);

% Image = crop_half_FOV(Image);

if tWeight ~= 0
    tNorm = abs(diff(Image,1,3));
    tNorm = tWeight * sum(tNorm(:));
    
    for i=1:size(bins,1)
        bin_temp = bins(i,:);
        Image_temp = Image(:,:,bin_temp,:,:);
        tNorm_temp = abs(diff(Image_temp,1,3));
        tNorm = tNorm + tWeight * sum(tNorm_temp(:));
    end
    
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

fNorm = fNorm/N;
tNorm = tNorm/N;
sNorm = sNorm/N;

Cost = sNorm + tNorm + fNorm;

if nargin == 3
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