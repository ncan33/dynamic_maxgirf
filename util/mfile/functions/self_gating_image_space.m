function bins = self_gating_image_space(Image)
Image = abs(Image);
[sx,sy,nof,set,sms] = size(Image);
signal = reshape(Image,[sx*sy,nof,set,sms]);
%signal(:,1:10,:,:) = [];% delete proton density images
signal = permute(signal,[1,4,2,3]);
signal = reshape(signal,[sx*sy*sms,nof,set]);

signal_mean = mean(signal,2);keyboard
signal_local_mean = (signal + circshift(signal,[0,1,0]) + circshift(signal,[0,-1,0]))/3;
signal_local_mean(:,1,:) = (signal(:,1,:) + signal(:,2,:))/2;
signal_local_mean(:,end,:) = (signal(:,end,:) + signal(:,end-1,:))/2;
signal_cardiac = abs(signal./signal_local_mean);

for i=1:set
    signal_mean_temp = signal_mean(:,:,i);
    signal_cardiac_temp = signal_cardiac(:,:,i);
    [~,order] = sort(signal_mean_temp,'descend');
    N = length(order);
    N = round(N/2);
    signal_cardiac_temp = signal_cardiac_temp(order(1:N),:);
    coeff = pca(signal_cardiac_temp.');
    cardiac_signal(:,i) = signal_cardiac_temp.'*coeff(:,1);

    cardiac_signal_temp = cardiac_signal(:,i);
    figure,plot(cardiac_signal_temp)
    
    temp = [cardiac_signal_temp,circshift(cardiac_signal_temp,1),circshift(cardiac_signal_temp,-1)];
    temp(1,2) = temp(1,1);
    temp(end,3) = temp(end,1);
    cardiac_signal_max = max(temp,[],2);
    cardiac_signal_min = min(temp,[],2);
    cardiac_signal_mid = (cardiac_signal_max+cardiac_signal_min)/2;
    temp = [cardiac_signal_max-cardiac_signal_temp,cardiac_signal_temp-cardiac_signal_min,abs(cardiac_signal_temp-cardiac_signal_mid)];
    [~,bins(:,i)] = min(temp,[],2);
end