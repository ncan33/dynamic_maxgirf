function disp = get_rigid_motion(Image,line,dim)
% Image: 3D image, x,y,t
% line:  which line to use
% dim:   the line is at whitch dim (1 or 2)
switch dim
    case 1
        line_profile = squeeze(Image(line,:,:));
    case 2
        line_profile = squeeze(Image(:,line,:));
end

diff_line_profile = diff(line_profile,1,1);
diff_sum = sum(abs(diff_line_profile),2);
diff_sum_smooth = conv(diff_sum,[1 1 1 1 1 1 1]);
[pks(:,1),pks(:,2)] = findpeaks(diff_sum_smooth);
pks = sortrows(pks,'descend');
heart_length = abs(pks(2,2) - pks(1,2))+1;

heart_kernel = ones(heart_length,1);
line_profile_conv = conv2(line_profile,heart_kernel);
[~,loc] = max(line_profile_conv,[],1);
disp = loc - round(mean(loc));