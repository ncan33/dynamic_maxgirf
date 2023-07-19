function Blind_Quant(Gd_curves,Mask_LV,Mask_MYO)



[sx,sy,sz,nof] = size(Gd_curves);
AIF = Gd_curves(repmat(Mask_LV,[1,1,1,nof]));
AIF = reshape(AIF,length(AIF)/nof,nof);
AIF = mean(AIF,1);
[~,peak_enhan_frame] = max(smooth(AIF));
[~,enhan_begin] = max(abs(diff(smooth(AIF))));
N_pre = round(enhan_begin*0.8);

Tissue_curves = Gd_curves(repmat(Mask_MYO,[1,1,1,nof]));
Tissue_curves = reshape(Tissue_curves,length(Tissue_curves)/nof,nof);

Nbins = 15;
bins = kmeans(Tissue_curves,Nbins);
for i=1:Nbins
    Tissue_curves_mean(i,:) = mean(Tissue_curves(bins==i,:));
end

idx_keep = true(1,Nbins);
idx_keep(min(Tissue_curves_mean(:,N_pre:end),[],2)<0) = false;
idx_keep(max(Tissue_curves_mean,[],2)==Tissue_curves_mean(:,peak_enhan_frame)) = false;
idx_keep(std(Tissue_curves_mean(:,1:N_pre),1,2)==0) = false;
idx_keep(std(Tissue_curves_mean(:,1:N_pre),1,2)>0.01) = false;

idx = false(size(bins));
for i=1:Nbins
    if idx_keep(i)
        idx(bins==i) = true;
    end
end

Tissue_curves = Tissue_curves(idx,:);

clear Tissue_curves_mean
Nbins = 10;
bins = kmeans(Tissue_curves,Nbins);
for i=1:Nbins
    Tissue_curves_mean(i,:) = mean(Tissue_curves(bins==i,:));
end


keyboard

AIF_para = Quant.fit_AIF(0:nof-1,double(AIF));
AIF_model = Quant.AIF_model(AIF_para,0:nof-1);


param = Fit_2CM(AIF_model,Tissue_curves_mean);
param = {param,0:nof-1};

for iter = 1:150
    iter
    options = optimoptions('lsqcurvefit','Display','off');
    AIF_update = lsqcurvefit(@Tissue_model,AIF_para,param,double(Tissue_curves_mean)',[],[],options);
    AIF_model = Quant.AIF_model(AIF_update,0:nof-1);
    param = Fit_2CM(AIF_model,Tissue_curves_mean);
    param = {param,0:nof-1};
    
    figure(2)
    clf
    plot(AIF_model)
    drawnow
end
param = param{1};