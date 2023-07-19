function w = density_comp_for_NUFFT(siz,NUFFT)
w = ones(siz);
for i=5
    tmp = NUFFT'* w;
    tmp = NUFFT * tmp;
    %w = 1./abs(tmp);
    %tmp_w = mean(mean(mean(abs(tmp),2),3),4);
    tmp_w = 1./real(tmp);
    %tmp_w = tmp_w/max(tmp_w(:));
    w = w.*tmp_w;
    %w = w.*repmat(tmp_w,[1 siz(2:end)]);
    %w = w/mean(w(:));
end
%w = w(:,1,1,1);