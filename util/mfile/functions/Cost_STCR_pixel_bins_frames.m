function [Cost_new,Cost,fNorm,tNorm,sNorm] = Cost_STCR_pixel_bins_frames(fUpdate, Image, Data, para, Cost_old)

%N = numel(Image);
sWeight = para.Recon.weight_sTV;
tWeight = para.Recon.weight_tTV;
bins = para.Recon.bins;

fNorm = sum(sum(sum(abs(fUpdate).^2,1),2),5);

nof = size(Image,3);

dt = abs(crop_half_FOV(zeros(size(Image),'like',Image)));
if tWeight ~= 0
    dt(:,:,1:nof-1,:,:) = abs(crop_half_FOV(Image(Data.Motion.idx_b) - Image(:,:,1:end-1,:,:)));
    dt(:,:,2:nof,:,:) = dt(:,:,2:nof,:,:) + abs(crop_half_FOV(Image(Data.Motion.idx_f) - Image(:,:,2:end,:,:)));
    dt(:,:,2:nof-1,:,:) = dt(:,:,2:nof-1,:,:)/2;
    tNorm = tWeight * sum(sum(sum(dt,1),2),5); clear dt

    for i=1:size(bins,1)
        bin_temp = bins(i,:);
        nof = sum(bin_temp);
        Image_temp = Image(:,:,bin_temp,:,:);
        Motion_temp = Data.Motion_bins{i};
        dt = abs(crop_half_FOV(zeros(size(Image_temp),'like',Image)));
        dt(:,:,1:nof-1,:,:) = abs(crop_half_FOV(Image_temp(Motion_temp.idx_b) - Image_temp(:,:,1:end-1,:,:)));
        dt(:,:,2:nof,:,:) = dt(:,:,2:nof,:,:) + abs(crop_half_FOV(Image_temp(Motion_temp.idx_f) - Image_temp(:,:,2:end,:,:)));
        dt(:,:,2:nof-1,:) = dt(:,:,2:nof-1,:)/2;
        
        tNorm(:,:,bins(i,:)) = tNorm(:,:,bins(i,:)) + tWeight * sum(sum(sum(dt,1),2),5); clear dt
    end
    
else
    tNorm = 0;
end

Image = crop_half_FOV(Image);

if sWeight ~= 0
    sx_norm = abs(diff(Image,1,2));
    sx_norm(:,end+1,:,:,:)=0;
    sy_norm = abs(diff(Image,1,1));
    sy_norm(end+1,:,:,:,:)=0;
    sNorm = sqrt(abs(sx_norm).^2+abs(sy_norm).^2);
    sNorm = sWeight * sum(sum(sum(sNorm,1),2),5);
else
    sNorm = 0;
end

Cost = sNorm + tNorm + fNorm;

if nargin == 4
    Cost_new = Cost;
    return
end

Cost_new = Cost_old;
fNorm = sum(fNorm);
tNorm = sum(tNorm);
sNorm = sum(sNorm);
Cost = sum(Cost);

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