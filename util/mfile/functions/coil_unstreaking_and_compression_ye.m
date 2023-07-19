function kSpace = coil_unstreaking_and_compression_ye(kSpace,kSpace_info)
idx = squeeze(sum(sum(sum(sum(kSpace)))) == 0);
kSpace(:,:,:,:,idx) = [];
kSpace = kSpace * 10^8;
kSpace_raw = kSpace;

kSpace = permute(kSpace,[1,2,4,5,3]);
[sx,nor,nof,ns,nc] = size(kSpace);
kSpace = reshape(kSpace,[sx,nor*nof,ns,nc]);

notwanted = kSpace_info.ray_mod==0;
kSpace(:,notwanted,:,:) = [];
phase_mod = kSpace_info.phase_mod(~notwanted);
angle_mod = kSpace_info.angle_mod(~notwanted);

nor_ref = length(angle_mod);
nof = nor_ref/nor;
nor_one = 50;
%ratio = nor_ref/nor_one;

%kSpace = kSpace(:,1:nor_ref,:,:);
kSpace = reshape(kSpace,[sx,nor_ref,1,ns,nc]);
phase_mod = phase_mod(1:nor_ref);
%angle_mod = angle_mod(1:nor_ref);

phase(1,:,1,1,1,1) = exp(-1i*phase_mod*0);
phase(1,:,1,1,1,2) = exp(-1i*phase_mod*2*pi/3);
phase(1,:,1,1,1,3) = exp(-1i*phase_mod*4*pi/3);

[kx,ky] = get_k_coor(sx,angle_mod.',0,round((sx+1)/2));
N = NUFFT.init_new(kx,ky,1.5,[4,4]);
im_ref = NUFFT.NUFFT_adj_new(kSpace,N);

kx = reshape(kx,[sx,nor,nof]);
ky = reshape(ky,[sx,nor,nof]);
kSpace = reshape(kSpace,[sx,nor,nof,ns,nc]);
N_one = NUFFT.init_new(kx,ky,1.5,[4,4]);
im_one = squeeze(NUFFT.NUFFT_adj_new(kSpace,N_one));

d = im_one-im_ref/nof;
nn = sum(sum(sum(abs(d).^2)));
dd = crop_half_FOV(im_ref);
dd = sum(sum(abs(dd).^2));

score = nn./dd;

for i=1:ns
    score_temp = squeeze(score(:,:,:,i,:));
    [~,order] = sort(score_temp,'descend');
    coils = order(1:8);
    d = im_one(:,:,:,i,coils)-im_ref(:,:,:,i,coils)/nof;
    k_d = NUFFT.NUFFT_new(d,N_one);
    H = 1-hann(288).^16;

end
keyboard
%for i=1:3
    %im_ref(:,:,:,:,i) = squeeze(NUFFT.NUFFT_adj_new(kSpace.*phase(:,:,:,:,:,i),N));


    %im_one(:,:,:,:,i) = squeeze(NUFFT.NUFFT_adj_new(kSpace(:,1:nor_one,:,:,:).*phase(:,1:nor_one,:,:,:,i),N_one));
    %im_one(:,:,:,:,i) = squeeze(NUFFT.NUFFT_adj_new(kSpace(:,1:nor_one,:,:,:),N_one));
%end
%max_one = max(max(abs(im_one)));
%max_ref = max(max(abs(im_ref)));




for i=1:ns
    for j=1:3
        im_one_temp = squeeze(im_one(:,:,i,:,j));
        im_ref_temp = squeeze(im_ref(:,:,i,:,j));
        scale_temp = max(max(abs(im_ref_temp)))./max(max(abs(im_one_temp)));
        n_temp = abs(im_one_temp.*scale_temp - im_ref_temp).^2;
        d_temp = abs(crop_half_FOV(im_ref_temp)).^2;
        score_temp(i,j,:) = sum(sum(n_temp))./sum(sum(d_temp));
    end
end

score = squeeze(sum(score_temp,2));
N = round(nc/2);

for i=1:3
    score_temp = score(i,:);
    score_temp = score_temp/min(score_temp);
    [~,order] = sort(score_temp);
    %threshold = score_temp(order==N);
    %if threshold<1.3
        threshold = 1.3;
    %end
    idx = score_temp>threshold;
    score(i,idx) = 1./score_temp(idx);
    score(i,~idx) = 1;
end

score = permute(score,[1,3,2]);

clear kSpace
for i=1:ns
    kSpace_temp = kSpace_raw(:,:,:,:,i).*score(i,1,:);
    kSpace_temp = permute(kSpace_temp,[1,2,4,3]);
    kSpace_temp = reshape(kSpace_temp,[sx*nor*nof,nc]);
    coeff = pca(kSpace_temp);
    kSpace_temp = kSpace_temp*coeff(:,1:8);
    kSpace_temp = reshape(kSpace_temp,[sx,nor,nof,8]);
    kSpace(:,:,:,:,i) = kSpace_temp;
end
