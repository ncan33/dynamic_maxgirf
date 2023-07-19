function plot_circle(position, radius, cmap)

x0 = position(1);
y0 = position(2);

theta = (1:360)/180*pi;
x = cos(theta) * radius + x0;
y = sin(theta) * radius + y0;

x1 = cos(theta) * radius/2 + x0;
y1 = sin(theta) * radius/2 + y0;

for i = 1:10:360
    plot_arrow_(x1(i), y1(i), x(i), y(i), cmap(i, :), 0.8, 2);
end
    
