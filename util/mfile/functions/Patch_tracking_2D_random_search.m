function llr = Patch_tracking_2D_random_search(Image,patch_size,search_size,patch_shift)

[sx,sy,nof] = size(Image);

if ~exist('patch_size')
    patch_size = [5,5];
end
if ~exist('search_size')
    search_size = [9,9];
end
if ~exist('patch_shift')
    patch_shift = [5,5];
end

x_begin_all = 1:patch_shift(1):sx-patch_size(1)+1;
y_begin_all = 1:patch_shift(2):sy-patch_size(2)+1;

if x_begin_all(end) + patch_size(1)-1 < sx
    x_begin_all(end+1) = sx-patch_size(1) + 1;
end
if y_begin_all(end) + patch_size(2)-1 < sy
    y_begin_all(end+1) = sy-patch_size(2) + 1;
end

Nx = length(x_begin_all);
Ny = length(y_begin_all);

std_map = std(Image,1,3);
int_map = sum(abs(Image),3);

offset_all = [];

N = 0;
for i=1:Nx
    for j=1:Ny
        
        x_temp = x_begin_all(i):x_begin_all(i)+patch_size(1)-1;
        y_temp = y_begin_all(j):y_begin_all(j)+patch_size(2)-1;
        
        [x_temp,y_temp] = ndgrid(x_temp,y_temp);
        idx_temp = sub2ind([sx,sy],x_temp,y_temp);
        
        N = N+1;
        %idx_all(N) = sum(vec(mask(x_temp,y_temp,z_temp)));
        idx_all(:,:,N) = idx_temp;
        int_patch(:,:,N) = int_map(idx_temp);
        std_patch(:,:,N) = std_map(idx_temp);
        
    end
end
    
    
idx = sum(sum(abs(int_patch)));
idx_std = sum(sum(abs(std_patch)));
keep = (idx_std > max(idx_std)/10) & (idx > max(idx)/10);
    
idx_all = idx_all(:,:,keep);
idx_add = (0:nof-1)*sx*sy;
idx_all = idx_all + permute(idx_add,[1,3,4,2]);
idx_all = idx_all(:,:,:);
Npatch = sum(keep(:));

[~,searching_order] = sort(rand([1,Npatch*nof]));

    
%% patch search
    
mask = zeros(size(Image));
n=0;

for i=1:Npatch*nof

    order_temp = searching_order(i);
    idx_temp = idx_all(:,:,order_temp);
    flag = sum(vec(mask(idx_temp) == 0));
        
    if flag
        n = n+1;
        [x_temp,y_temp,c_temp] = ind2sub([sx,sy,nof],idx_temp);
        x_begin = x_temp(1);
        y_begin = y_temp(1);
        c_temp = c_temp(1);
           
        mask(idx_temp) = mask(idx_temp) + 1;
        offset_all(:,c_temp) = [x_begin,y_begin];
        if c_temp~=1
            for icycle = c_temp-1:-1:1
                
                x_window_begin = max(x_begin-round((search_size(1)-patch_size(1))/2),1); x_window_end = min(x_window_begin-1+search_size(1),sx);
                if x_window_begin==1
                    x_window_end = max(x_window_end,x_window_begin+patch_size(1)-1);
                end
                if x_window_end==sx
                    x_window_begin = min(x_window_begin,x_window_end-patch_size(1)+1);
                end
                
                y_window_begin = max(y_begin-round((search_size(2)-patch_size(2))/2),1); y_window_end = min(y_window_begin-1+search_size(2),sy);
                if y_window_begin==1
                    y_window_end = max(y_window_end,y_window_begin+patch_size(2)-1);
                end
                if y_window_end==sy
                    y_window_begin = min(y_window_begin,y_window_end-patch_size(2)+1);
                end
                
                patch_temp = Image(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,icycle+1);
                search_window = Image(x_window_begin:x_window_end,y_window_begin:y_window_end,icycle);
                
                [offset] = patch_search_one_2D(patch_temp,search_window);
                
                begin_all(1) = x_window_begin-1+offset(1); x_begin = begin_all(1);
                begin_all(2) = y_window_begin-1+offset(2); y_begin = begin_all(2);
                
                offset_all(:,icycle) = begin_all;
                mask(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,icycle) = mask(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,icycle) + 1;
            end
        end
        if c_temp~=nof
            x_begin = x_temp(1);
            y_begin = y_temp(1);
            
            for icycle = c_temp:1:nof-1
                
                x_window_begin = max(x_begin-round((search_size(1)-patch_size(1))/2),1); x_window_end = min(x_window_begin-1+search_size(1),sx);
                if x_window_begin==1
                    x_window_end = max(x_window_end,x_window_begin+patch_size(1)-1);
                end
                if x_window_end==sx
                    x_window_begin = min(x_window_begin,x_window_end-patch_size(1)+1);
                end
                
                y_window_begin = max(y_begin-round((search_size(2)-patch_size(2))/2),1); y_window_end = min(y_window_begin-1+search_size(2),sy);
                if y_window_begin==1
                    y_window_end = max(y_window_end,y_window_begin+patch_size(2)-1);
                end
                if y_window_end==sy
                    y_window_begin = min(y_window_begin,y_window_end-patch_size(2)+1);
                end
                
                patch_temp = Image(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,icycle);
                search_window = Image(x_window_begin:x_window_end,y_window_begin:y_window_end,icycle+1);
                
                [offset] = patch_search_one_2D(patch_temp,search_window);
                
                begin_all(1) = x_window_begin-1+offset(1); x_begin = begin_all(1);
                begin_all(2) = y_window_begin-1+offset(2); y_begin = begin_all(2);
                
                offset_all(:,icycle+1) = begin_all;
                mask(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,icycle+1) = mask(x_begin:x_begin+patch_size(1)-1,y_begin:y_begin+patch_size(2)-1,icycle+1) + 1;
            end
        end

        offset_all_phase(:,:,n) = offset_all;
    end
end
%%
% n = sum(phase_all);

patch_begin_all = offset_all_phase;
patch_end_all = offset_all_phase + reshape(patch_size(1:2),2,1) - 1;

idx = zeros([patch_size,nof,n]);
for i=1:n
    %     fprintf([num2str(i/n),'\n'])
    for j=1:nof
        x_temp = patch_begin_all(1,j,i):patch_end_all(1,j,i);
        y_temp = patch_begin_all(2,j,i):patch_end_all(2,j,i);
        
        [x_temp,y_temp,p_temp] = ndgrid(x_temp,y_temp,j);
        idx_temp = sub2ind([sx,sy,nof],x_temp,y_temp,p_temp);
        idx(:,:,j,i) = idx_temp;
    end
end
idx = reshape(idx,prod(patch_size),nof,n);

mask_LR = mask;
mask_LR(mask_LR==0) = 1;
mask_LR = 1./mask_LR;
mask = logical(mask);
Npatch = n;

llr.idx = idx;
llr.Npatch = Npatch;
llr.mask_intensity = mask_LR;
llr.mask = mask;

return










