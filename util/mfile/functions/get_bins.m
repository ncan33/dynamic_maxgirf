function bin = get_bins(cardiac_bins,respiration_bins)
card_N = max(cardiac_bins);
resp_N = max(respiration_bins);

N = card_N*resp_N;

card_n = repmat(1:card_N,[1,resp_N]);
resp_n = repmat(1:resp_N,[1,card_N]);

for i=1:N
    bin(i,:) = cardiac_bins == card_n(i) & respiration_bins == resp_n(i);
end