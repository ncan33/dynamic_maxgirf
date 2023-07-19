function llr = Patch_tracking_3D_with_guide(Data,para,mask,offset_init,patch_size,search_size,patch_shift)
if isfield(Data,'field')
    im_bin = Data.first_guess;
else
    im_bin = Data;
end
[sx,sy,nof,nslice] = size(im_bin);

Nphase = size(para.Recon.bins,1);
Ncycle = nof/Nphase;
im_bin = reshape(im_bin,[sx,sy,Nphase,Ncycle,nslice]);
im_bin = permute(im_bin,[1,2,5,3,4]);

if ~exist('patch_size')
    patch_size = [5,5,2];
end
if ~exist('search_size')
    search_size = [9,9,4];
end
if ~exist('patch_shift')
    patch_shift = [3,3,1];
end

% Nx = floor(sx/patch_size(1));
% Ny = floor(sx/patch_size(2));
% Nz = ceil(nslice/patch_size(3));
% Nc = ceil(Nphase/patch_size(4));

x_begin_all = 1:patch_shift(1):sx-patch_size(1)+1;
y_begin_all = 1:patch_shift(2):sy-patch_size(2)+1;
z_begin_all = 1:patch_shift(3):nslice-patch_size(3)+1;
% c_begin_all = 1:Nphase:Nphase-patch_size(4)+1;

if z_begin_all(end) + patch_size(3)-1 < nslice
    z_begin_all(end+1) = nslice-patch_size(3) + 1;
end
if x_begin_all(end) + patch_size(1)-1 < sx
    x_begin_all(end+1) = sx-patch_size(1) + 1;
end
if y_begin_all(end) + patch_size(2)-1 < sy
    y_begin_all(end+1) = sy-patch_size(2) + 1;
end

Nx = length(x_begin_all);
Ny = length(y_begin_all);
Nz = length(z_begin_all);
% Nc = length(c_begin_all);

%Npatch_all = Nx*Ny*Nz*Nc;
%int_patch = zeros([patch_size,Npatch_all]);
%std_map = std(im_bin,1,5);
%int_map = sum(abs(im_bin),5);
%std_patch = zeros([patch_size,Npatch_all]);
N = 0;


for i=1:Nx
    for j=1:Ny
        for k=1:Nz
%             for l=1:Nc
                x_temp = x_begin_all(i):x_begin_all(i)+patch_size(1)-1;
                y_temp = y_begin_all(j):y_begin_all(j)+patch_size(2)-1;
                z_temp = z_begin_all(k):z_begin_all(k)+patch_size(3)-1;
%                 c_temp = c_begin_all(l):c_begin_all(l)+patch_size(4)-1;
                N = N+1;
                idx_all(N) = sum(vec(mask(x_temp,y_temp,z_temp)));
                %int_patch(:,:,:,:,N) = int_map(x_temp,y_temp,z_temp,c_temp);
                %std_patch(:,:,:,:,N) = std_map(x_temp,y_temp,z_temp,c_temp);
%             end
        end
    end
end

%idx = sum(sum(sum(sum(abs(int_patch)))));
%idx_std = sum(sum(sum(sum(abs(std_patch)))));
%keep = (idx_std > max(idx_std)/10) & (idx > max(idx)/10);
keep = idx_all~=0;

%% patch search

%offset_init = squeeze(sum(offset_init,2));
offset_init = reshape(offset_init,[3,Nphase,Ncycle]);
offset_init = offset_init - offset_init(:,:,1);
offset_init(:,:,2:end) = diff(offset_init,1,3);
% offset_init = squeeze(mean(offset_init,2));
% offset_init(3,:) = offset_init(3,:)/3;
offset_init = round(offset_init);
% offset_init = cat(1,offset_init,zeros(1,Ncycle));
% offset_init = permute(offset_init,[2,1]);
offset_init = -offset_init;

Npatch = sum(keep(:));
%patch_series_all = zeros([patch_size,Ncycle,Npatch]);
offset_all = zeros(Ncycle,Nphase,3,Npatch);
% N = 0;

% idx_keep = find(keep);
% [x,y,z,c] = ind2sub([Nx,Ny,Nz,Nc],idx_keep);
x_all = zeros([patch_size,Nphase,Npatch]);
y_all = zeros([patch_size,Nphase,Npatch]);
z_all = zeros([patch_size,Nphase,Npatch]);
% c_all = zeros([patch_size,Npatch]);
% for i=1:Npatch
%     tic
%     fprintf([num2str(i),'/',num2str(Npatch),' '])
%     x_temp = x_begin_all(x(i)):x_begin_all(x(i))+patch_size(1)-1;
%     y_temp = y_begin_all(y(i)):y_begin_all(y(i))+patch_size(2)-1;
%     z_temp = z_begin_all(z(i)):z_begin_all(z(i))+patch_size(3)-1;
%     c_temp = c_begin_all(c(i)):c_begin_all(c(i))+patch_size(4)-1;
%     [x_all(:,:,:,:,i),y_all(:,:,:,:,i),z_all(:,:,:,:,i),c_all(:,:,:,:,i)] = ndgrid(x_temp,y_temp,z_temp,c_temp);
%     [offset,~] = patch_search_4D(im_bin,[x_temp(1),y_temp(1),z_temp(1),c_temp(1)],patch_size,patch_size+4);
%     offset_all(:,:,i) = offset;
%     toc
% end




for nphase = 1:Nphase
    N = 0;
    im_temp = squeeze(im_bin(:,:,:,nphase,:));
    offset_temp = squeeze(offset_init(:,nphase,:));
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
%             for l=1:Nc
                N = N+1;
                if keep(N)
                   
%                     tic
                    n = find(find(keep)==N);
%                     fprintf([num2str(n),'/',num2str(Npatch),' '])
                    x_temp = x_begin_all(i):x_begin_all(i)+patch_size(1)-1;
                    y_temp = y_begin_all(j):y_begin_all(j)+patch_size(2)-1;
                    z_temp = z_begin_all(k):z_begin_all(k)+patch_size(3)-1;
%                     c_temp = c_begin_all(l):c_begin_all(l)+patch_size(4)-1;
                    
                    [offset,~] = patch_search_3D_with_guide(im_temp,[x_temp(1),y_temp(1),z_temp(1)],patch_size,search_size,offset_temp);
                    %[offset,patch_series]=patch_search_4D(im_bin,[x_temp(1),y_temp(1),z_temp(1),c_temp(1)],patch_size,patch_size+4);
                    %patch_series_all(:,:,:,:,:,n) = patch_series;
                    offset_all(:,nphase,:,n) = offset;
                    [x_all(:,:,:,nphase,n),y_all(:,:,:,nphase,n),z_all(:,:,:,nphase,n)] = ndgrid(x_temp,y_temp,z_temp);
%                     toc
                end
%             end
        end
    end
end
end

% idx_drop = squeeze(sum(sum(abs(offset_all),1),3)) == 0;
% offset_all(:,:,idx_drop) = [];
% x_all(:,:,:,:,idx_drop) = [];
% y_all(:,:,:,:,idx_drop) = [];
% z_all(:,:,:,:,idx_drop) = [];
% c_all(:,:,:,:,idx_drop) = [];
% Npatch = size(offset_all,3);


for i=Ncycle:-1:1
    offset_all(i,:,:) = sum(offset_all(1:i,:,:),1);
end
offset_all = permute(offset_all,[3,5,6,2,1,4]);
x_all = permute(x_all,[1,2,3,4,6,5]) + offset_all(1,:,:,:,:,:);
y_all = permute(y_all,[1,2,3,4,6,5]) + offset_all(2,:,:,:,:,:);
z_all = permute(z_all,[1,2,3,4,6,5]) + offset_all(3,:,:,:,:,:);
% c_all = permute(c_all,[1,2,3,4,6,5]) + offset_all(4,:,:,:,:,:);
idx = sub2ind([sx,sy,nslice,Nphase],x_all,y_all,z_all);

idx_add = reshape((0:Ncycle*Nphase-1)*sx*sy*nslice,[1,Ncycle*Nphase]);

% idx = idx+idx_add;
idx = reshape(idx,prod(patch_size),Nphase*Ncycle,Npatch);
idx = idx+idx_add;

% patch_series_all = im_bin(idx);
% for i=1:Npatch
%     [u,s,v] = svd(patch_series_all(:,:,i),0);
%     s = s-s(end);
%     s(s<0) = 0;
%     patch_all_LR(:,:,i) = u*s*v';
% end

% im_LR = zeros(sx*sx*nslice*Nphase,Ncycle);
mask = zeros(sx,sx,nslice,Nphase*Ncycle,'single');
for i=1:Npatch
    % im_LR(idx(:,:,i)) = im_LR(idx(:,:,i)) + patch_all_LR(:,:,i);
    mask(idx(:,:,i)) = mask(idx(:,:,i)) + 1;
end

mask_LR = mask;
mask_LR(mask_LR==0) = 1;
mask_LR = 1./mask_LR;
mask_LR = reshape(mask_LR,[sx,sy,nslice,Nphase*Ncycle]);
% im_LR = im_LR./mask_LR;

mask = logical(mask);
mask = reshape(mask,[sx,sy,nslice,Nphase*Ncycle]);

llr.idx = idx;
llr.Npatch = Npatch;
llr.mask_intensity = mask_LR;
llr.mask = mask;

end


