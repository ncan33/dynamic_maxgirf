function Data = off_resonance_low_rank_approx(Data, para)

[sx, sy, nof] = size(Data.off_res.map);

f_image = zeros([sx, sy, nof, para.L], 'single');
f_kspace = zeros([length(para.t), 1, nof, para.L], 'single');

for iframe = 1:nof
    fmap = Data.off_res.map(:, :, iframe);
    E_mat = exp(-1i * 2 * pi * para.t' * fmap(:)');
    [count, center] = hist(fmap(:), para.n_hist_bins);
    ek = exp(-1i * 2 * pi * para.t' * center);
    [U, ~, ~] = svd(repmat(sqrt(count), [size(ek,1) 1]) .* ek, 'econ');
    
    B_mat = U(:, 1:para.L);
    C_mat = pinv(B_mat) * E_mat;
    C_mat = reshape(C_mat, [para.L, sx, sy]);
   
    f_image(:, :, iframe, :) = permute(C_mat,[2, 3, 1]);
    f_kspace(:, :, iframe, :) = B_mat;
end

Data.off_res.f_im = f_image;
Data.off_res.f_k = f_kspace;