function update = patch_low_rank_only(Image,llr)

Image = permute(Image,[1,2,4,3]);
patch_all = Image(llr.idx);
update = zeros(size(Image),'like',Image);

% update_tv = compute_tTV_yt(permute(patch_all,[1,3,2]),para.Recon.weight_tTV,1e-7);
% update_tv = permute(update_tv,[1,3,2]);
% update_t = zeros(size(Image),'like',Image);

for i=1:llr.Npatch
%     fprintf([num2str(i/llr.Npatch),'\n'])
    [u,s,v] = svd(patch_all(:,:,i),0);
    s = s-s(end);
    s(s<0) = 0;
    update(llr.idx(:,:,i)) = update(llr.idx(:,:,i)) + u*s*v';
%     update_t(llr.idx(:,:,i)) = update_t(llr.idx(:,:,i)) + update_tv(:,:,i);
    %patch_all(:,:,i) = u*s*v'; 
end

update = update.*llr.mask_intensity;
% update = update_t.*llr.mask_intensity + update*0.1;    
update = permute(update,[1,2,4,3]);