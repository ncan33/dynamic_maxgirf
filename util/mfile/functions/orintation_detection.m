% 
% UCAIR, University of Utah
% Ye Tian, phye1988@gmail.com
%
% Orin = orintation_detection(image)
%
% Get orintation flag for input image.
%
%
% Input:
%      Image - should be scalar, 2D image.

function Orin = orintation_detection(image)

o1 = image;
o2 = rot90(image,1);
o3 = rot90(image,2);
o4 = rot90(image,3);
o5 = flipud(o1);
o6 = flipud(o2);
o7 = flipud(o3);
o8 = flipud(o4);
    
fig1 = figure;
subplot(2,4,1); imagesc(o1); colormap gray; brighten(0.4); axis image; axis off; title '1';
subplot(2,4,2); imagesc(o2); colormap gray; brighten(0.4); axis image; axis off; title '2';
subplot(2,4,3); imagesc(o3); colormap gray; brighten(0.4); axis image; axis off; title '3';
subplot(2,4,4); imagesc(o4); colormap gray; brighten(0.4); axis image; axis off; title '4';
subplot(2,4,5); imagesc(o5); colormap gray; brighten(0.4); axis image; axis off; title '5';
subplot(2,4,6); imagesc(o6); colormap gray; brighten(0.4); axis image; axis off; title '6';
subplot(2,4,7); imagesc(o7); colormap gray; brighten(0.4); axis image; axis off; title '7';
subplot(2,4,8); imagesc(o8); colormap gray; brighten(0.4); axis image; axis off; title '8';

prompt = 'Which one is the correct orintation? ';
Orin = input(prompt);
close(fig1);