function [Line,fig] = line_profile(Image,line,dim)

if ~isreal(Image)
    Image = abs(Image);
end

switch dim
    case 1
        Image = squeeze(Image(line,:,:));
    case 2
        Image = squeeze(Image(:,line,:));
end

if size(Image,1) == 1 || size(Image,2) == 1
    fig = figure;
    plot(Image)
else
    fig = figure;
    imagesc(Image)
    axis image
    colormap gray
end

Line = Image;

return
