function play_mri_video(n_frames, fps, video_matrix)
    % Note: fps is equal to the acquistion fps (2 fps) divided by the
    % number of interleaves for the full kspace trajectory (10) times the
    % number of arms per frame (5), which should be undersampled. In this
    % example calculation, the fps for 5 arms per frame should be 2/10*5 = 
    % 4 fps.
    
    % set n_frames to 'all' if you want to play every frame
    
    if isnumeric(n_frames == 'all')
        n_frames = size(video_matrix, 3);
    end
    
    for i = 1:n_frames
        a = abs(video_matrix(:,:,i));
        a = fliplr(rot90(a, -1));
        imagesc(a); axis image; colorbar; colormap gray
        caxis([0 abs(mean(squeeze(max(squeeze(max(video_matrix))))))])
        
        pause(1/fps)
    end
end