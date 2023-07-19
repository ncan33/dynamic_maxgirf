function [bins,signal_all] = self_gating_image_space_2(Image)
Image = abs(Image);
[sx,sy,nof,set,sms] = size(Image);
signal = reshape(Image,[sx*sy,nof,set,sms]);
%signal(:,1:10,:,:) = [];% delete proton density images
signal = permute(signal,[1,4,2,3]);
signal = reshape(signal,[sx*sy*sms,nof,set]);

signal_mean = mean(signal,2);
signal_cardiac_temp = signal(:,:,1);
sum_signal = sum(signal_cardiac_temp,1);
[~,idx] = max(sum_signal);
for i=1:set
    signal_mean_temp = signal_mean(:,:,i);
    signal_cardiac_temp = signal(:,:,i);
    [~,order] = sort(signal_mean_temp,'descend');
    N = length(order);
    N = round(N/2);
    signal_cardiac_temp = signal_cardiac_temp(order(1:N),:);
    sum_signal = sum(signal_cardiac_temp,1);
    %[~,idx] = max(sum_signal);
    [~,order] = sort(signal_cardiac_temp(:,idx),'descend');
    N = round(N/10*2);
    signal_cardiac_temp = signal_cardiac_temp(order(1:N),:);
    
    signal_std = std(signal_cardiac_temp(:,idx-5:idx+5).');
    [~,order] = sort(signal_std,'descend');
    N = round(N/4);
    signal_cardiac_temp = signal_cardiac_temp(order(1:N),:);

    for j=1:N
        temp = signal_cardiac_temp(j,:);
        [mins,locs] = findpeaks(-temp);
        mins = -mins;
        cali = interp1(locs,mins,(1:nof));
        cali(1:locs(1)) = min(temp(1:locs(1)),cali(locs(1)));
        cali(locs(end):end) = min(temp(locs(end):end),cali(locs(end)));
        cali_temp = temp - cali;
        cali_temp(cali_temp>0) = 0;
        while sum(cali_temp)<0
            [~,locs_add] = findpeaks(-cali_temp);
            locs = [locs,locs_add];
            mins = [mins,temp(locs_add)];
            [locs,order] = sort(locs);
            mins = mins(order);
            cali = interp1(locs,mins,(1:nof));
            cali(1:locs(1)) = min(temp(1:locs(1)),cali(locs(1)));
            cali(locs(end):end) = min(temp(locs(end):end),cali(locs(end)));
            cali_temp = temp - cali;
            cali_temp(cali_temp>0) = 0;
        end
        temp = temp - cali;
        
        [maxs,locs] = findpeaks(temp);
        cali = interp1(locs,maxs,(1:nof));
        cali(1:locs(1)) = max(temp(1:locs(1)),cali(locs(1)));
        cali(locs(end):end) = max(temp(locs(end):end),cali(locs(end)));
        cali_temp = cali - temp;
        cali_temp(cali_temp>0) = 0;
        while sum(cali_temp)<0
            [~,locs_add] = findpeaks(-cali_temp);
            locs = [locs,locs_add];
            maxs = [maxs,temp(locs_add)];
            [locs,order] = sort(locs);
            maxs = maxs(order);
            cali = interp1(locs,maxs,(1:nof));
            cali(1:locs(1)) = max(temp(1:locs(1)),cali(locs(1)));
            cali(locs(end):end) = max(temp(locs(end):end),cali(locs(end)));
            cali_temp = cali - temp;
            cali_temp(cali_temp>0) = 0;
        end
        signal_cardiac_temp(j,:) = temp./cali;
    end

    %signal_all(:,i) = sum(signal_cardiac_temp,1);
    coeff = pca(signal_cardiac_temp.');
    signal_all(:,i) = signal_cardiac_temp.' * coeff(:,1);
end

signal_all = signal_all - min(signal_all,[],1);

for i=1:set
    temp = signal_all(:,i);
    temp_smoothed = smooth(temp);
    [mins,locmins] = findpeaks(-temp);
    mins = -mins;
    idx_notwanted = mins > temp_smoothed(locmins);
    mins(idx_notwanted) = [];
    locmins(idx_notwanted) = [];
    min_temp = interp1(locmins,mins,(1:nof).');
    min_temp(1:locmins(1)) = min(temp(1:locmins(1)),min_temp(locmins(1)));
    min_temp(locmins(end):end) = min(temp(locmins(end):end),min_temp(locmins(end)));
    
    flag_temp = temp - min_temp;
    flag_temp(flag_temp>0) = 0;
    
    while sum(flag_temp) < 0
        [~,locs_add] = findpeaks(-flag_temp);
        locmins = [locmins;locs_add];
        mins = [mins;temp(locs_add)];
        [locmins,order] = sort(locmins);
        mins = mins(order);
        min_temp = interp1(locmins,mins,(1:nof).');
        min_temp(1:locmins(1)) = min(temp(1:locmins(1)),min_temp(locmins(1)));
        min_temp(locmins(end):end) = min(temp(locmins(end):end),min_temp(locmins(end)));
        flag_temp = temp - min_temp;
        flag_temp(flag_temp>0) = 0;
    end
    
    [maxs,locmaxs] = findpeaks(temp);
    idx_notwanted = maxs < temp_smoothed(locmaxs);
    maxs(idx_notwanted) = [];
    locmaxs(idx_notwanted) = [];
    max_temp = interp1(locmaxs,maxs,(1:nof).');
    max_temp(1:locmaxs(1)) = max(temp(1:locmaxs(1)),max_temp(locmaxs(1)));
    max_temp(locmaxs(end):end) = max(temp(locmaxs(end):end),max_temp(locmaxs(end)));

    flag_temp = max_temp - temp;
    flag_temp(flag_temp>0) = 0;
    
    while sum(flag_temp) < 0
        [~,locs_add] = findpeaks(-flag_temp);
        locmaxs = [locmaxs;locs_add];
        maxs = [maxs;temp(locs_add)];
        [locmaxs,order] = sort(locmaxs);
        maxs = maxs(order);
        max_temp = interp1(locmaxs,maxs,(1:nof).');
        max_temp(1:locmaxs(1)) = max(temp(1:locmaxs(1)),max_temp(locmaxs(1)));
        max_temp(locmaxs(end):end) = max(temp(locmaxs(end):end),max_temp(locmaxs(end)));
        flag_temp = max_temp - temp;
        flag_temp(flag_temp>0) = 0;
    end
    
    mid_temp = (max_temp+min_temp)/2;
    
    signal_all(:,i) = (temp-min_temp)./(max_temp-min_temp);
    
    d = [temp - min_temp, temp - max_temp, temp - mid_temp];

    [~,bins(:,i)] = min(abs(d),[],2);
end


