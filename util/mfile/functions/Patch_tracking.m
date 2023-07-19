function Data = Patch_tracking(Data,para)

% [sx,sy,nof,~,nset,nSMS] = size(Data.first_guess);
% Nphase = size(para.Recon.bins,1);
% Ncycle = nof/Nphase;
% im_bin = reshape(Data.first_guess,[sx,sy,Nphase,Ncycle,nSMS*nset]);
% 
% im_bin = permute(im_bin,[1,2,5,3,4]);
% nslice = nSMS*nset;
% order = 1:nSMS:nslice;
% for i=nSMS:-1:2
%     order = [order,i:nset:nslice];
% end
% [~,order_back] = sort(order);
% im_bin = im_bin(:,:,order,:,:);
% 
% [x,y,z,ph,cy] = ndgrid(1:sx,1:sy,1:nslice,1:Nphase,1:Ncycle);
% [xv,yv,zv,phv,cyv] = ndgrid(1:sx,1:sy,1:(nslice-1)/(nslice*7/1.8-1):nslice,1:Nphase,1:Ncycle);
% im_bin = interpn(x,y,z,ph,cy,im_bin,xv,yv,zv,phv,cyv);
% 
% nslice = size(im_bin,3);
[sx,sy,nof,nslice] = size(Data.first_guess);
Nphase = size(para.Recon.bins,1);
Ncycle = nof/Nphase;
im_bin = reshape(Data.first_guess,[sx,sy,Nphase,Ncycle,nslice]);
im_bin = permute(im_bin,[1,2,5,3,4]);

patch_size = [7,7,2,Nphase];
search_size = [11,11,4,Nphase];
patch_shift = [4,4,1];


% Nx = floor(sx/patch_size(1));
% Ny = floor(sx/patch_size(2));
% Nz = ceil(nslice/patch_size(3));
% Nc = ceil(Nphase/patch_size(4)); 

x_begin_all = 1:patch_shift(1):sx-patch_size(1)+1;
y_begin_all = 1:patch_shift(2):sy-patch_size(2)+1;
z_begin_all = 1:patch_shift(3):nslice-patch_size(3)+1;
c_begin_all = 1:Nphase:Nphase-patch_size(4)+1;

if z_begin_all(end) + patch_size(3)-1 < nslice
    z_begin_all(end+1) = nslice-patch_size(3) + 1;
end
keyboard
Nx = length(x_begin_all);
Ny = length(y_begin_all);
Nz = length(z_begin_all);
Nc = length(c_begin_all);

Npatch_all = Nx*Ny*Nz*Nc;
int_patch = zeros([patch_size,Npatch_all]);
std_map = std(im_bin,1,5);
int_map = sum(abs(im_bin),5);
std_patch = zeros([patch_size,Npatch_all]);
N = 0;
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            for l=1:Nc
                x_temp = x_begin_all(i):x_begin_all(i)+patch_size(1)-1;
                y_temp = y_begin_all(j):y_begin_all(j)+patch_size(2)-1;
                z_temp = z_begin_all(k):z_begin_all(k)+patch_size(3)-1;
                c_temp = c_begin_all(l):c_begin_all(l)+patch_size(4)-1;
                N = N+1;
                int_patch(:,:,:,:,N) = int_map(x_temp,y_temp,z_temp,c_temp);
                std_patch(:,:,:,:,N) = std_map(x_temp,y_temp,z_temp,c_temp);
            end
        end
    end
end

idx = sum(sum(sum(sum(abs(int_patch)))));
idx_std = sum(sum(sum(sum(abs(std_patch)))));
keep = (idx_std > max(idx_std)/10) & (idx > max(idx)/10);


%% patch search
Npatch = sum(keep(:));
%patch_series_all = zeros([patch_size,Ncycle,Npatch]);
offset_all = zeros(Ncycle,4,Npatch);
% N = 0;

% idx_keep = find(keep);
% [x,y,z,c] = ind2sub([Nx,Ny,Nz,Nc],idx_keep);
x_all = zeros([patch_size,Npatch]);
y_all = zeros([patch_size,Npatch]);
z_all = zeros([patch_size,Npatch]);
c_all = zeros([patch_size,Npatch]);
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

N = 0;

for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            for l=1:Nc
                N = N+1;
                if keep(N)
                    tic
                    n = find(find(keep)==N);
                    fprintf([num2str(n),'/',num2str(Npatch),' '])
                    x_temp = x_begin_all(i):x_begin_all(i)+patch_size(1)-1;
                    y_temp = y_begin_all(j):y_begin_all(j)+patch_size(2)-1;
                    z_temp = z_begin_all(k):z_begin_all(k)+patch_size(3)-1;
                    c_temp = c_begin_all(l):c_begin_all(l)+patch_size(4)-1;
                    
                    [offset,~]=patch_search_4D(im_bin,[x_temp(1),y_temp(1),z_temp(1),c_temp(1)],patch_size,search_size);
                    %[offset,patch_series]=patch_search_4D(im_bin,[x_temp(1),y_temp(1),z_temp(1),c_temp(1)],patch_size,patch_size+4);
                    %patch_series_all(:,:,:,:,:,n) = patch_series;
                    offset_all(:,:,n) = offset;
                    [x_all(:,:,:,:,n),y_all(:,:,:,:,n),z_all(:,:,:,:,n),c_all(:,:,:,:,n)] = ndgrid(x_temp,y_temp,z_temp,c_temp);
                    toc
                end
            end
        end
    end
end


idx_drop = squeeze(sum(sum(abs(offset_all)))) == 0;
offset_all(:,:,idx_drop) = [];
x_all(:,:,:,:,idx_drop) = [];
y_all(:,:,:,:,idx_drop) = [];
z_all(:,:,:,:,idx_drop) = [];
c_all(:,:,:,:,idx_drop) = [];
Npatch = size(offset_all,3);



for i=Ncycle:-1:1
    offset_all(i,:,:) = sum(offset_all(1:i,:,:),1);
end
offset_all = permute(offset_all,[2,4,5,6,1,3]);
x_all = permute(x_all,[1,2,3,4,6,5]) + offset_all(1,:,:,:,:,:);
y_all = permute(y_all,[1,2,3,4,6,5]) + offset_all(2,:,:,:,:,:);
z_all = permute(z_all,[1,2,3,4,6,5]) + offset_all(3,:,:,:,:,:);
c_all = permute(c_all,[1,2,3,4,6,5]) + offset_all(4,:,:,:,:,:);
idx = sub2ind([sx,sy,nslice,Nphase],x_all,y_all,z_all,c_all);

idx_add = reshape((0:Ncycle-1)*sx*sy*nslice*Nphase,[1,1,1,1,Ncycle]);
idx = idx+idx_add;
idx = reshape(idx,prod(patch_size),Ncycle,Npatch);


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

Data.llr.idx = idx;
Data.llr.Npatch = Npatch;
Data.llr.mask_intensity = mask_LR;
Data.llr.mask = mask;


