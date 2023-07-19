function show_compare_yt(im1,im2)
figure
nof = size(im1,3);
i = 1;
while i
    j = mod(i,nof);
    if j == 0
        j = nof;
    end
    subplot(2,1,1)
    imagesc(abs(im1(:,:,j)))
    colormap gray
    axis image
    brighten(0.3)
    title(num2str(j))
    drawnow
    
    subplot(2,1,2)
    imagesc(abs(im2(:,:,j)))
    colormap gray
    axis image
    brighten(0.3)
    i = i+1;
    pause
end