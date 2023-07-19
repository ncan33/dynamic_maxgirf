function llr = Patch_init_with_guide(Data,para,mask,offset_init,patch_size,patch_shift)

[sx,sy,nof,nslice] = size(Data.first_guess);
Nphase = size(para.Recon.bins,1);
Ncycle = nof/Nphase;
%im_bin = reshape(Data.first_guess,[sx,sy,Nphase,Ncycle,nslice]);
%im_bin = permute(im_bin,[1,2,5,3,4]);

if ~exist('patch_size')
    patch_size = [5,5,5];
end
%search_size = [5,5,5,Nphase];
if ~exist('patch_shift')
    patch_shift = 3;
end


% Nx = floor(sx/patch_size(1));
% Ny = floor(sx/patch_size(2));
% Nz = ceil(nslice/patch_size(3));
% Nc = ceil(Nphase/patch_size(4));

x_begin_all = 1:patch_shift:sx-patch_size(1)+1;
y_begin_all = 1:patch_shift:sy-patch_size(2)+1;
z_begin_all = 1:patch_shift:nslice-patch_size(3)+1;
%c_begin_all = 1:Nphase:Nphase-patch_size(4)+1;

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
%Nc = length(c_begin_all);

%Npatch_all = Nx*Ny*Nz*Nc;
%int_patch = zeros([patch_size,Npatch_all]);
%std_map = std(im_bin,1,5);
%int_map = sum(abs(im_bin),5);
%std_patch = zeros([patch_size,Npatch_all]);
N = 0;


for i=1:Nx
    for j=1:Ny
        for k=1:Nz
%            for l=1:Nc
                x_temp = x_begin_all(i):x_begin_all(i)+patch_size(1)-1;
                y_temp = y_begin_all(j):y_begin_all(j)+patch_size(2)-1;
                z_temp = z_begin_all(k):z_begin_all(k)+patch_size(3)-1;
                %c_temp = c_begin_all(l):c_begin_all(l)+patch_size(4)-1;
                N = N+1;
                idx_all(N) = sum(vec(mask(x_temp,y_temp,z_temp)));
                %int_patch(:,:,:,:,N) = int_map(x_temp,y_temp,z_temp,c_temp);
                %std_patch(:,:,:,:,N) = std_map(x_temp,y_temp,z_temp,c_temp);
%            end
        end
    end
end

%idx = sum(sum(sum(sum(abs(int_patch)))));
%idx_std = sum(sum(sum(sum(abs(std_patch)))));
%keep = (idx_std > max(idx_std)/10) & (idx > max(idx)/10);
keep = idx_all~=0;

%% patch search

%offset_init = squeeze(sum(offset_init,2));
%offset_init = reshape(offset_init,[3,Nphase,Ncycle]);
offset_init = offset_init - offset_init(:,1);
%offset_init(:,2:end) = diff(offset_init,1,2);
%offset_init = squeeze(mean(offset_init,2));
%offset_init(3,:) = offset_init(3,:)/3;
offset_init = round(offset_init);
offset_init(:,2:end) = diff(offset_init,1,2);
%offset_init = cat(1,offset_init,zeros(1,Nphase,Ncycle));
offset_init = permute(offset_init,[2,1]);
offset_init = -offset_init;

Npatch = sum(keep(:));
%patch_series_all = zeros([patch_size,Ncycle,Npatch]);
offset_all = zeros(nof,3,Npatch);
% N = 0;

% idx_keep = find(keep);
% [x,y,z,c] = ind2sub([Nx,Ny,Nz,Nc],idx_keep);
x_all = zeros([patch_size,Npatch]);
y_all = zeros([patch_size,Npatch]);
z_all = zeros([patch_size,Npatch]);
%c_all = zeros([patch_size,Npatch]);
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
                    
                    %[offset,~]=patch_search_4D_with_guide(im_bin,[x_temp(1),y_temp(1),z_temp(1),c_temp(1)],patch_size,search_size,offset_init);
                    offset = offset_init;
                    
                    %[offset,patch_series]=patch_search_4D(im_bin,[x_temp(1),y_temp(1),z_temp(1),c_temp(1)],patch_size,patch_size+4);
                    %patch_series_all(:,:,:,:,:,n) = patch_series;
                    offset_all(:,:,n) = offset;
                    [x_all(:,:,:,n),y_all(:,:,:,n),z_all(:,:,:,n)] = ndgrid(x_temp,y_temp,z_temp);
%                     toc
                end
%             end
        end
    end
end

% idx_drop = squeeze(sum(sum(abs(offset_all)))) == 0;
% offset_all(:,:,idx_drop) = [];
% x_all(:,:,:,:,idx_drop) = [];
% y_all(:,:,:,:,idx_drop) = [];
% z_all(:,:,:,:,idx_drop) = [];
% c_all(:,:,:,:,idx_drop) = [];
Npatch = size(offset_all,3);

for i=nof:-1:1
    offset_all(i,:,:) = sum(offset_all(1:i,:,:),1);
end

offset_all = permute(offset_all,[2,4,5,1,3]);
x_all = permute(x_all,[1,2,3,5,4]) + offset_all(1,:,:,:,:,:);
y_all = permute(y_all,[1,2,3,5,4]) + offset_all(2,:,:,:,:,:);
z_all = permute(z_all,[1,2,3,5,4]) + offset_all(3,:,:,:,:,:);

z_all(z_all<1) = 1;
z_all(z_all>nslice) = nslice;
%c_all = permute(c_all,[1,2,3,4,6,5]) + offset_all(4,:,:,:,:,:);
idx = sub2ind([sx,sy,nslice],squeeze(x_all),squeeze(y_all),squeeze(z_all));

idx_add = reshape((0:nof-1)*sx*sy*nslice,[1,1,1,nof]);
idx = idx+idx_add;
idx = reshape(idx,prod(patch_size),nof,Npatch);


% patch_series_all = im_bin(idx);
% for i=1:Npatch
%     [u,s,v] = svd(patch_series_all(:,:,i),0);
%     s = s-s(end);
%     s(s<0) = 0;
%     patch_all_LR(:,:,i) = u*s*v';
% end

% im_LR = zeros(sx*sx*nslice*Nphase,Ncycle);
mask = zeros(sx,sx,nslice,nof,'single');
for i=1:Npatch
    % im_LR(idx(:,:,i)) = im_LR(idx(:,:,i)) + patch_all_LR(:,:,i);
    mask(idx(:,:,i)) = mask(idx(:,:,i)) + 1;
end
mask_LR = mask;
mask_LR(mask_LR==0) = 1;
mask_LR = 1./mask_LR;
mask_LR = reshape(mask_LR,[sx,sy,nslice,nof]);
% im_LR = im_LR./mask_LR;

mask = logical(mask);
mask = reshape(mask,[sx,sy,nslice,nof]);

llr.idx = idx;
llr.Npatch = Npatch;
llr.mask_intensity = mask_LR;
llr.mask = mask;

end



function [offset,patch_series] = patch_search_4D_with_guide(image,patch_loc,patch_size,search_box,offset_init)

siz = size(image);

patch_s = patch_loc;
patch_e = patch_loc + patch_size - 1;

% patch_s_0 = patch_s;
% patch_e_0 = patch_e;

offset = real(zeros(siz(5),4,'like',image));
%offset = offset_init;

patch_series = zeros([patch_size,siz(5)],'like',image);

for nof = 1:siz(5)-1
    patch_s = patch_s + offset(nof,:);
    patch_e = patch_e + offset(nof,:);
    
    patch_temp = image(patch_s(1):patch_e(1),patch_s(2):patch_e(2),patch_s(3):patch_e(3),patch_s(4):patch_e(4),nof);
    patch_series(:,:,:,:,nof) = patch_temp;
    
    %%%
    %patch_temp = patch_temp - mean(patch_temp(:));
    %%%
    
    search_box_s = round(patch_s + (patch_size-1)/2 - (search_box-1)/2) + offset_init(nof+1,:);
    search_box_e = round(patch_s + (patch_size-1)/2 + (search_box-1)/2) + offset_init(nof+1,:);
    
    search_box_s = max([1,1,1,1],search_box_s);
    search_box_e = min(siz(1:4),search_box_e);
    
    search_box_temp = image(search_box_s(1):search_box_e(1),search_box_s(2):search_box_e(2),search_box_s(3):search_box_e(3),search_box_s(4):search_box_e(4),nof+1);
    search_box_size = size(search_box_temp);

    Npatches = search_box_size - patch_size + 1;
    
    rela_loc = patch_s - search_box_s + 1;
    
    patch_all = zeros([patch_size,prod(Npatches)],'like',image);
    %distance = zeros(1,prod(Npatches),'like',image);
    for s1 = 1:Npatches(1)
        for s2 = 1:Npatches(2)
            for s3 = 1:Npatches(3)
                for s4 = 1:Npatches(4)
                    n = sub2ind(Npatches,s1,s2,s3,s4);
                    patch_all(:,:,:,:,n) = search_box_temp(s1:s1+patch_size(1)-1,s2:s2+patch_size(2)-1,s3:s3+patch_size(3)-1,s4:s4+patch_size(4)-1);
                    %distance(n) = sqrt(sum(([s1,s2,s3,s4] - rela_loc).^2));
                end
            end
        end
    end
    
    %%%
    %patch_all = patch_all - mean(mean(mean(mean(patch_all))));
    %%%
    
    d = abs(patch_temp - patch_all);
    d = reshape(d,prod(patch_size),prod(Npatches));
    d = sum(d);
    %d = d + distance*sum(patch_temp(:))/100;
    [~,idx] = min(d);
    [offset(nof+1,1),offset(nof+1,2),offset(nof+1,3),offset(nof+1,4)] = ind2sub(Npatches,idx);
    offset(nof+1,:) = offset(nof+1,:) - rela_loc;
end
nof = nof+1;
patch_s = patch_s + offset(nof,:);
patch_e = patch_e + offset(nof,:);
patch_temp = image(patch_s(1):patch_e(1),patch_s(2):patch_e(2),patch_s(3):patch_e(3),patch_s(4):patch_e(4),nof);
patch_series(:,:,:,:,nof) = patch_temp;

end


