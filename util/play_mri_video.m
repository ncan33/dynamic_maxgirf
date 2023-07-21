function play_mri_video(n_frames, fps, video_matrix, save_video)
    % This function is called as follows:
    %
    % play_mri_video(n_frames, fps, video_matrix, save_video)
    % 
    % - n_frames: number of frame to show. set n_frames to 'all' if you
    %   want to play every frame.
    % - fps: frames per second
    % - video_matrix: 3D matrix where the 3rd dimension is the temporal one
    % - save_video: set this to 1 if you wish to save video in the
    %   make_video directory's tmp folder for videomaking
    %
    % Note: fps is equal to the acquistion fps (2 fps) divided by the
    % number of interleaves for the full kspace trajectory (10) times the
    % number of arms per frame (5), which should be undersampled. In this
    % example calculation, the fps for 5 arms per frame should be 2/10*5 = 
    % 4 fps.
    
    arguments
        n_frames
        fps
        video_matrix
        save_video = 0
    end
    
    if ~isnumeric(n_frames)
        n_frames = size(video_matrix, 3);
    end
    
    if save_video
        figure;
        input(['Resize the figure as you wish. When done, press enter', ...
            ' in command window'])
        
        vidName = 'mri_video';
        video = VideoWriter(vidName);
        video.FrameRate = fps;
        open(video)
    end
    
    for i = 1:n_frames
        %a = abs(video_matrix(:,:,i));
        a = abs(video_matrix(70:349,  70:349, i)); % for NUFFT
        a = fliplr(rot90(a, -1)); % for NUFFT
        imagesc(a); axis image; colorbar; colormap gray
        caxis([0 abs(mean(squeeze(max(squeeze(max(video_matrix))))))])
        
        %% play vs. save video
        if ~save_video
            pause(1/fps)
        else
            saveas(gcf, ['/server/home/ncan/make_video/tmp', num2str(i), '.png'])
            image = imread(['/server/home/ncan/make_video/tmp', num2str(i), '.png']);
            writeVideo(video, image);
            disp(['Frame ', num2str(i),' done!'])
        end
    end
    
    if save_video
        disp('Frames successfully saved in ncan/make_video_tmp')
        convert_command = sprintf('ffmpeg -y -i %s.avi %s.mp4', vidName, vidName);
        system(convert_command)
    end
end