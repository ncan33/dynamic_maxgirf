function slice = pick_slice(Image)

Image = squeeze(abs(Image));
%Image = permute(Image,[1,3,2]);
siz = size(Image);
Image = reshape(Image,[siz(1),siz(2)*siz(3)]);

temp = figure;
imagesc(Image)
axis image
axis off
colormap gray
brighten(0.4)

prompt = 'Which slice you want to extract cardiac signal from? ';
slice = input(prompt);
close(temp)