function normalized_signal = normalize_to_0_1(signal)
signal = signal(:);
N = length(signal);
signal = signal - min(signal);

[mins,locs] = findpeaks(-signal);
mins = -mins;
cali = interp1(locs,mins,(1:N).');
cali(1:locs(1)) = min(signal(1:locs(1)),cali(locs(1)));
cali(locs(end):end) = min(signal(locs(end):end),cali(locs(end)));
cali_temp = signal - cali;
cali_temp(cali_temp>0) = 0;
while sum(cali_temp)<0
    [~,locs_add] = findpeaks(-cali_temp);
    locs = [locs;locs_add];
    mins = [mins;signal(locs_add)];
    [locs,order] = sort(locs);
    mins = mins(order);
    cali = interp1(locs,mins,(1:N).');
    cali(1:locs(1)) = min(signal(1:locs(1)),cali(locs(1)));
    cali(locs(end):end) = min(signal(locs(end):end),cali(locs(end)));
    cali_temp = signal - cali;
    cali_temp(cali_temp>0) = 0;
end

cali_min = cali;
        
[maxs,locs] = findpeaks(signal);
cali = interp1(locs,maxs,(1:N).');
cali(1:locs(1)) = max(signal(1:locs(1)),cali(locs(1)));
cali(locs(end):end) = max(signal(locs(end):end),cali(locs(end)));
cali_temp = cali - signal;
cali_temp(cali_temp>0) = 0;
while sum(cali_temp)<0
    [~,locs_add] = findpeaks(-cali_temp);
    locs = [locs;locs_add];
    maxs = [maxs;signal(locs_add)];
    [locs,order] = sort(locs);
    maxs = maxs(order);
    cali = interp1(locs,maxs,(1:N).');
    cali(1:locs(1)) = max(signal(1:locs(1)),cali(locs(1)));
    cali(locs(end):end) = max(signal(locs(end):end),cali(locs(end)));
    cali_temp = cali - signal;
    cali_temp(cali_temp>0) = 0;
end

normalized_signal = (signal-cali_min)./(cali-cali_min);