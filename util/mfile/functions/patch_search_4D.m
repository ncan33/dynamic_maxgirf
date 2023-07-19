function [offset,patch_series] = patch_search_4D(image,patch_loc,patch_size,search_box)

siz = size(image);

patch_s = patch_loc;
patch_e = patch_loc + patch_size - 1;

% patch_s_0 = patch_s;
% patch_e_0 = patch_e;

offset = real(zeros(siz(5),4,'like',image));
patch_series = zeros([patch_size,siz(5)],'like',image);
for nof = 1:siz(5)-1
    patch_s = patch_s + offset(nof,:);
    patch_e = patch_e + offset(nof,:);
    
    patch_temp = image(patch_s(1):patch_e(1),patch_s(2):patch_e(2),patch_s(3):patch_e(3),patch_s(4):patch_e(4),nof);
    patch_series(:,:,:,:,nof) = patch_temp;
    
    %%%
    %patch_temp = patch_temp - mean(patch_temp(:));
    %%%
    
    search_box_s = round(patch_loc + (patch_size-1)/2 - (search_box-1)/2);
    search_box_e = round(patch_loc + (patch_size-1)/2 + (search_box-1)/2);
    
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
