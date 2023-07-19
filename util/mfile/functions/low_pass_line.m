function kSpace_lp = low_pass_line(kSpace,idx)

siz = size(idx);
kSpace_line = kSpace(idx);

% [kSpace_line(1,:),lp_filt] = lowpass(kSpace_line(1,:),0.42);
% lp_filt = double(lp_filt);
for i=1:siz(1)
    kSpace_line(i,:) = lowpass(kSpace_line(i,:),0.25);
%     kSpace_line(i,:) = fftfilt(lp_filt,double(kSpace_line(i,:)));
end

kSpace_lp = zeros(size(kSpace));
mask = kSpace_lp;
for i=1:siz(1)
    for j=1:siz(2)
        kSpace_lp(idx(i,j)) = kSpace_lp(idx(i,j)) + kSpace_line(i,j);
        mask(idx(i,j)) = mask(idx(i,j)) + 1;
    end
end

imagesc(mask==0)
mask(mask==0) = 1;
mask = 1./mask;
kSpace_lp = kSpace_lp.*mask;

idx_zero = abs(kSpace_lp)==0;
kSpace_lp(idx_zero) = kSpace(idx_zero);