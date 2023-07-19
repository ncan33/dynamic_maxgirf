function [xcoor, ycoor] = get_k_coor_RING(sx,theta,ifNUFFT,kCenter,S)

xcoor = (1:sx) - kCenter;
ycoor = xcoor;
xcoor = bsxfun(@times,xcoor',cos(theta));
ycoor = bsxfun(@times,ycoor',sin(theta));

if ifNUFFT == 1
    xcoor = xcoor/sx;
    ycoor = ycoor/sx;
end

correction = S.*permute([cos(theta);sin(theta)],[4,1,2,3]);
correction = squeeze(sum(correction,2));
xcoor = xcoor + correction(1,:,:);
ycoor = ycoor + correction(1,:,:);

end
