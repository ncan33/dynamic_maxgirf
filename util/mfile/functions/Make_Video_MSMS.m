function Make_Video_MSMS(file_name)

file = dir(file_name);
load(file.name)
Image = reshape_MSMS(Image);

nof = size(Image,3);

video_name = [file.name(1:end-4),'.avi'];

v = VideoWriter(video_name);
v.FrameRate = 10;
v.Quality = 100;
open(v)

Video = figure;

for i=1:nof
    imagesc(Image(:,:,i))
    colormap gray
    axis image
    axis off
    brighten(0.4)
    drawnow
    image = getframe;
    writeVideo(v,image)

end

close(v)
close(Video)
