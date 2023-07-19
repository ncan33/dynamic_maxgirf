function image = sliding_window_recon(Data, nor_one, nor_sl)

[sx, nor, nof, nc] = size(Data.kSpace);
nsms = size(Data.phase_mod, 5);
nor_total = nor*nof;
nof_sl = floor((nor_total - nor_one) / nor_sl) + 1;
if nof_sl == 0
    nof_sl = 1;
end

kSpace = reshape(Data.kSpace, [sx, nor*nof, 1, nc]);
kx = reshape(Data.kx, [sx, nor*nof]);
ky = reshape(Data.ky, [sx, nor*nof]);
phase_mod = reshape(Data.phase_mod, [1, nor*nof, 1, 1, nsms]);

image = zeros([sx, sx, nof_sl, nc, nsms], 'like', Data.kSpace);

for i=1:nof_sl
    fprintf(sprintf('%g/%g time frame\n', i, nof_sl))
    idx = (1:nor_one) + nor_sl * (i-1);
    
    N = NUFFT.init(kx(:, idx), ky(:, idx), 1, [6, 6], sx, sx);
    image(:, :, i, :, :) = NUFFT.NUFFT_adj(kSpace(:, idx, :, :) .* conj(phase_mod(:, idx, :, :, :)), N);
end