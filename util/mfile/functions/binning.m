function bins = binning(signal,nbins)

signal_sorted = sort(signal);
n = length(signal)/nbins;
threshold(1) = 0;
bins = zeros(size(signal));
threshold = zeros(1,nbins);
for i=1:nbins-1
    threshold(i+1) = signal_sorted(round(n*i));
    bins(signal<threshold(i+1) & signal>=threshold(i)) = i;
end
bins(bins==0) = nbins;
