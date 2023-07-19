function w = density_comp_for_GNUFFT(siz,G)
w = ones(siz);
for i=5
    tmp = GROG.GNUFFT_adj(w,G);
    tmp = GROG.GNUFFT(tmp,G);
    %w = 1./abs(tmp);
    tmp_w = mean(mean(mean(abs(tmp),2),3),4);
    tmp_w = 1./tmp_w;
    tmp_w = tmp_w/max(tmp_w(:));
    
    w = w.*repmat(tmp_w,[1 siz(2:end)]);
    %w = w/mean(w(:));
end
w = w(:,1,1,1);