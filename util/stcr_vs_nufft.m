function stcr_vs_nufft(frame, im_echo_1, im_echo_2, NUFFT_im_echo_1, NUFFT_im_echo_2)
    % Pass 'all' or 'end' to the parameter frame. In the future, the
    % ability to inspect a specific frame will be introduced.
    if frame == 'end'
        figure;
        subplot(3,1,1); imagesc(abs([im_echo_1(:,:,end), im_echo_2(:,:,end)])); axis image;
        colormap gray; colorbar
        title('Echo 1 (left) and Echo 2 (right) of the final frame of STCR')
        NUFFT_1 = NUFFT_im_echo_1(:,:,end);
        NUFFT_2 = NUFFT_im_echo_2(:,:,end);

        sz = size(NUFFT_1);
        sz_stcr = size(im_echo_1);

        lo = sz(1)/2 - sz_stcr(1)/2;
        hi = sz(1)/2 + sz_stcr(1)/2 - 1;

        NUFFT_1 = NUFFT_1(lo:hi, lo:hi);
        NUFFT_2 = NUFFT_2(lo:hi, lo:hi);
        
        subplot(3,1,2); imagesc(abs([fliplr(rot90(NUFFT_1,-1)), fliplr(rot90(NUFFT_2,-1))])); axis image;
        colormap gray; colorbar
        title('Echo 1 (left) and Echo 2 (right) of the final frame of NUFFT')
        
    elseif frame == 'all'
        error('Function is under development. This comparison can only be made for the last frame.')
    end
end