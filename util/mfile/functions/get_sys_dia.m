function [idx_s,idx_d] = get_sys_dia(Image,mask)
% [sys_loc,dia_loc] = get_sys_dia(Image)

Image = abs(Image);

if ~exist('mask')
    temp = figure;
    imagesc(sum(Image,3));
    colormap gray
    brighten(0.4)
    axis image
    axis off
    h = imfreehand(gca,'closed',false); 
    mask = createMask(h);
    close(temp)
end


cardiac_signal = squeeze(sum(sum(Image.*mask)));

[~,loc_s] = findpeaks(-cardiac_signal);
[~,loc_d] = findpeaks(cardiac_signal);

nof = size(Image,3);
idx = true(1,nof);

idx_d = false(1,nof);
idx_s = false(1,nof);

idx_d(loc_d) = true;
idx_s(loc_s) = true;

idx(loc_s) = false;
idx(loc_d) = false;

f = smooth((1:nof)',cardiac_signal,0.03,'loess');

for i=1:nof
    if idx(i)
        if cardiac_signal(i) > f(i)
            idx_d(i) = true;
        else
            idx_s(i) = true;
        end
    end
end