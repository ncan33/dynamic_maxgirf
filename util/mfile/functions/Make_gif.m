function Make_gif(file_name, Image)
Nframe = size(Image,3);
h = figure;
% Image = Image./max(vec(Image))*255;
for n = 1:Nframe
    imagesc(Image(:,:,n))
    axis image
    colormap gray
    brighten(0.2)
    axis off
    drawnow
    set(gcf,'Position',[100,100,400,400])
    set(gca, 'pos', [0,0,1,1]);
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
% imind = (Image(:,:,n));
    % Write to the GIF File
    if n == 1
        imwrite(imind,cm,file_name,'gif','DelayTime',0.2, 'Loopcount',inf);
    else
        imwrite(imind,cm,file_name,'gif','DelayTime',0.2,'WriteMode','append');
    end
end