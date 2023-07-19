function textMRI(grid, Text, positionx, positiony, offsetx, offsety)
if nargin == 2
    positionx = 'center';
    positiony = 'top';
    offsetx = 0;
    offsety = 0;
elseif nargin == 3
    positiony = 'top';
    offsetx = 0;
    offsety = 0;
elseif nargin == 4
    offsetx = 0;
    offsety = 0;
elseif nargin == 5
    offsety = 0;
end


x = (1:grid(1)) * 1/grid(1);
switch positionx
    case 'center'
        x = x - 1/grid(1)/2;
    case 'left'
        x = x - 1/grid(1);
    case 'right'
        
end
x = x + offsetx;

y = (grid(2):-1:1) * 1/grid(2);
switch positiony
    case 'top'
       
    case 'bottom'
        y = y - 1/grid(2);
   
end
y = y + offsety;

[x, y] = meshgrid(x, y);
x = vec(x');
y = vec(y');

for i = 1:length(Text)
    text(x(i), y(i), Text{i}, 'Color', 'white', 'units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 20)
end