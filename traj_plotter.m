clf
hold on
axis square
xlim([-0.5 0.5])
ylim([-0.5 0.5])

for j = 1:20
    for i = 0:62
        plot(kx((10*i+1):(10*i+1)+9,j), ky((10*i+1):(10*i+1)+9,j));
        title(['Interleaf: ', num2str(j), '. Sample: ', num2str(i*10)]);
        pause(0.001);
    end
end

%{
for i = 1:100
    plot(kx(:,i), ky(:,i));
    title(num2str(i));
    pause(0.2);
end
%}

%{
for j = 1:10
    for i = 0:89
        plot(kx((100*i+1):(100*i+1)+99,j), ky((100*i+1):(100*i+1)+99,j));
        title(['Interleaf: ', num2str(j), '. Sample: ', num2str(i*100)]);
        pause(0.01);
    end
end
%}