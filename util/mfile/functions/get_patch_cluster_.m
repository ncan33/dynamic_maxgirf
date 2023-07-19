function ttv = get_patch_cluster_(Image, patch_size, shift, n_cluster)
tic
nof = size(Image, 3);


fun = @(x) reshape(kmeans(reshape(x.data, [9, nof]).', n_cluster), [1, 1, nof]);
B = blockproc(abs(Image), [3, 3], fun);

fun2 = @(x) patch_cluster_tv(x, B);
toc

tic
ttv = blockproc(Image, [3, 3], fun2);
toc