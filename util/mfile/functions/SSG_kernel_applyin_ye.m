function SB = SSG_kernel_applyin_ye(krnl, MB, options)

h = options.patch_height;
w = options.patch_width;

[ny,nx,nc,~] = size(MB);
%ref_sms = sum(ref,4);

blks_MB = zeros(h*w,(ny-h+1)*(nx-w+1),nc);

for c=1:nc
    blks_MB(:,:,c) = im2col(MB(:,:,c,1),[h w],'sliding');
end

blks_MB = reshape(blks_MB,h,w,ny-h+1,nx-w+1,nc);
blks_MB = permute(blks_MB,[3 4 1 2 5]);

S = blks_MB;
S = reshape(S,size(S,1)*size(S,2),[]);

% resolved = multiprod(S,krnl,[1 2]);
% resolved = permute(resolved,[2 1 3 4 5]);
resolved = S * squeeze(krnl);
resolved = permute(resolved, [3,1,2]);
resolved = reshape(resolved,[],nx-w+1,size(resolved,3),size(resolved,4),size(resolved,5));
recon = zeros(size(resolved,1)+h-1,size(resolved,2)+w-1,size(resolved,3),size(resolved,4),size(resolved,5));
recon(floor((h+1)/2):size(resolved,1)+floor((h+1)/2-1 ) ,...
    floor((w+1)/2):size(resolved,2)+floor((w+1)/2-1 ),:,:,:)=resolved;
% SB = zeros(size(recon,1),size(recon,2),size(recon,3),size(recon,4));

SB = recon;
%% start test
% figure,imagesc(sos(fftshift2(ifft2(SB)))), axis image, colormap gray
% 
% K = reshape(krnl,[6,6,8,8]);
% K_large = zeros(288,288,8,8);
% K_large(142:147, 142:147, :,:) = K;
% K = K_large;
% K = fftshift2(K);
% K = ifft2(K,288,288);
% K = fftshift2(K);
% 
% o3 = rot90(K,2);
% 
% 
% MB_test = fftshift2(ifft2(fftshift2(MB)));
% 
% temp = squeeze(sum(MB_test.*o3,3));
% figure,imagesc(sos(temp)), axis image, colormap gray



end