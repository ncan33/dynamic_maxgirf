function para = auto_gating_cardiac_SMS_yt(Image_ref,para)
%{
short_long_flag = isfield(para{1}.Recon,'long_axis_idx');

for i=1:size(Image_ref,4)
    if short_long_flag && i == para{i}.Recon.long_axis_idx
        sys = auto_gating_long_axis(Image_ref(:,:,:,i,2));
    else
        sys = auto_gating_short_axis(Image_ref(:,:,:,i,2));
    end

    para{i}.Recon.bins(1,:) = sys;
    para{i}.Recon.bins(2,:) = ~sys;
end
%}

cardiac_signal = self_gating_image_space_4(crop_half_FOV(Image_ref));

for i=1:size(Image_ref,4)
    para{i}.Recon.bins(1,:) = cardiac_signal(:,i)==1;
    para{i}.Recon.bins(2,:) = cardiac_signal(:,i)==2;
    para{i}.Recon.bins(3,:) = cardiac_signal(:,i)==3;
    
    idx = sum(para{i}.Recon.bins,2)<=2;
    para{i}.Recon.bins(idx,:) = [];
    
    para{i}.cardiac_signal = cardiac_signal;
end
