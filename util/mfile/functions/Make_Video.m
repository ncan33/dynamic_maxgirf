function Make_Video(video_name, image)

%file = dir(file_name);
%load(file.name)
%Image = reshape_MSMS(Image);

nof = size(image,3);

%video_name = [file.name(1:end-4),'.avi'];
%{
v = VideoWriter(file_name, 'Uncompressed AVI');
% v.FrameRate = 1000/5.58/2; % frames per second
% v.FrameRate = 1000/5.56/10;
% v.FrameRate = 1000/6.002/2;
% v.FrameRate = 10;
v.FrameRate = 1000 / 40.11 / 1;

% v.Quality = 50;
open(v)

Video = figure(1);

for i=1:nof
    imagesc(Image(:,:,i))
    colormap gray
    axis image
    axis off
%     brighten(0.2)
    drawnow
    set(gca,'pos',[0.05 0.05 0.9 0.9])
    image = getframe;
    writeVideo(v,image)

end

close(v)
close(Video)
%}

v = vision.VideoFileWriter(video_name, 'FileFormat', 'avi');
    v.FrameRate = 1000 / 40.11 / 1;
    
    for ii = 1:nof
        v(image(:,:,ii) * 2.5);
    end
    clear im_all
    release(v)
    
    convert_command = sprintf('ffmpeg -y -i %s.avi -pix_fmt yuv420p %s.mp4', video_name(1:end-4), video_name(1:end-4));
    
    system(convert_command)
