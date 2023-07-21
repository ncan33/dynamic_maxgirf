function play_mri_video(n_frames, fps, video_matrix, save_video, resize_figure, ...
    convert_to_mov, vid_name)
    % This function is called as follows:
    %
    % play_mri_video(n_frames, fps, video_matrix, save_video, resize figure)
    % 
    % - n_frames: number of frame to show. set n_frames to 'all' if you
    %   want to play every frame.
    % - fps: frames per second
    % - video_matrix: 3D matrix where the 3rd dimension is the temporal one
    % - save_video: set this to 1 if you wish to save video in the
    %   make_video directory's tmp folder for videomaking
    % - resize_figure: if set to 1, a figure and the user is asked to
    %   resize it before plotting. if set to 0, the video will be plotted
    %   in the current figure (if none exists, a standard size figure will
    %   be created)
    %
    % Note: fps is equal to the acquistion fps (2 fps) time the number of
    % interleaves for the full kspace trajectory (10) divided by the number
    % of arms per frame (5), which should be undersampled. In this example
    % calculation, the fps for 5 arms per frame should be 2*10/5 = 4 fps.
    
    arguments
        n_frames
        fps
        video_matrix
        save_video = 0
        resize_figure = 0
        convert_to_mov = 1
        vid_name = 'mri_video'
    end
    
    if ~isnumeric(n_frames)
        n_frames = size(video_matrix, 3);
    end
    
    if save_video
        if resize_figure
            figure;
            input(['Resize the figure as you wish. When done, press ', ...
                'enter in command window'])
        end
        
        video = VideoWriter(vid_name);
        video.FrameRate = fps;
        open(video)
    end
    
    for i = 1:n_frames
        a = abs(video_matrix(:,:,i));
        %a = abs(video_matrix(70:349,  70:349, i)); % for NUFFT
        %a = fliplr(rot90(a, -1)); % for NUFFT
        imagesc(a); axis image; colorbar; colormap gray
        caxis([0 abs(mean(squeeze(max(squeeze(max(video_matrix))))))])
        %caxis([0 0.5])
        
        %% play vs. save video
        if ~save_video
            pause(1/fps)
        else
            saveas(gcf, ['/server/home/ncan/make_video/tmp/', num2str(i), '.png'])
            image = imread(['/server/home/ncan/make_video/tmp/', num2str(i), '.png']);
            writeVideo(video, image);
            disp(['Frame ', num2str(i),' done!'])
        end
    end
    
    if save_video && convert_to_mov
        close(video)
        disp('Frames successfully saved in ncan/make_video_tmp')
        convert_command = sprintf('ffmpeg -i %s.avi -c:v copy -c:a copy %s.mov', vid_name, vid_name);
        system(convert_command)
        system(['rm ',vid_name,'.avi'])
    end
end