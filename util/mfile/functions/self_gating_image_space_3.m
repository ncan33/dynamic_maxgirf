function [bins,signal_all] = self_gating_image_space_3(Image)
Image = abs(Image);
[sx,sy,nof,set,sms] = size(Image);
signal = reshape(Image,[sx*sy,nof,set,sms]);
signal = permute(signal,[1,4,2,3]);
signal = reshape(signal,[sx*sy*sms,nof,set]);
signal_max = max(signal,[],2);
N = size(signal_max,1);
N = round(N*0.2);
for i=1:set
    [~,order] = sort(signal_max(:,:,i),'descend');
    cardiac_signal(:,:,i) = signal(order(1:N),:,i);
end

signal_max = max(cardiac_signal,[],2);
cardiac_signal = cardiac_signal./signal_max;
pre_contrast_max = max(cardiac_signal(:,1:20,:),[],2);
for i=1:set
    threshold = 0.1;
    cardiac_signal_temp = cardiac_signal(:,:,i);
    temp = cardiac_signal_temp(pre_contrast_max(:,:,i)<threshold,:);
    while size(temp,1) > round(N*0.25)
        threshold = threshold * 0.95;
        temp = cardiac_signal_temp(pre_contrast_max(:,:,i)<threshold,:);
    end
    while size(temp,1) < round(N*0.025)
        threshold = threshold * 1.2;
        temp = cardiac_signal_temp(pre_contrast_max(:,:,i)<threshold,:);
    end
    cardiac_signal_temp = temp.';
    for j=1:size(cardiac_signal_temp,2)
        cardiac_signal_temp(:,j) = normalize_to_0_1(cardiac_signal_temp(:,j));
    end
    coeff = pca(cardiac_signal_temp);
    signal_temp = cardiac_signal_temp*coeff(:,1);

    signal_all(:,i) = normalize_to_0_1_2(signal_temp);
end

bins = zeros(size(signal_all));
bins(signal_all<0.2) = 1;
bins(signal_all>0.8) = 2;
bins(bins==0) = 3;
