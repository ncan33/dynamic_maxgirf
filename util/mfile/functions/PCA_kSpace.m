function kSpace = PCA_kSpace(kSpace)

siz = size(kSpace);
ND = length(siz);
kSpace = double(reshape(kSpace,prod(siz(1:ND-1)),siz(ND)));
scale_kspace = 10^8;
kSpace = kSpace*scale_kspace;
coeff = pca(kSpace);
if size(coeff,2)>8
    kSpace = kSpace*coeff(:,1:8);
    kSpace = reshape(kSpace,[siz(1:ND-1),8]);
else
    kSpace = reshape(kSpace,[siz(1:ND-1),size(coeff,2)]);
end
kSpace = permute(kSpace,[1 2 3 5 4]);
kSpace = single(kSpace);

end