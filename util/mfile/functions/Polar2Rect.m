function RectImage = Polar2Rect(Image,sx,sy,X,Y)

[x,y] = meshgrid(1:sx,1:sy);
warning off
RectImage = griddata(X,Y,double(Image),x,y,'nearest');