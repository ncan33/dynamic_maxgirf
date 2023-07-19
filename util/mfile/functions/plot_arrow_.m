function [h1, h2, h3] = plot_arrow_( x1,y1,x2,y2,colorin, alphain, linewidth )
%
% plot_arrow - plots an arrow to the current plot
%
% format:   handles = plot_arrow( x1,y1,x2,y2 [,options...] )
%
% input:    x1,y1   - starting point
%           x2,y2   - end point
%           options - come as pairs of "property","value" as defined for "line" and "patch"
%                     controls, see matlab help for listing of these properties.
%                     note that not all properties where added, one might add them at the end of this file.
%                     
%                     additional options are:
%                     'headwidth':  relative to complete arrow size, default value is 0.07
%                     'headheight': relative to complete arrow size, default value is 0.15
%                     (encoded are maximal values if pixels, for the case that the arrow is very long)
%
% output:   handles - handles of the graphical elements building the arrow
%
% Example:  plot_arrow( -1,-1,15,12,'linewidth',2,'color',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5] );
%           plot_arrow( 0,0,5,4,'linewidth',2,'headwidth',0.25,'headheight',0.33 );
%           plot_arrow;   % will launch demo


% =============================================
% constants (can be edited)
% =============================================
alpha_      = 0.3;   % head length
beta        = 0.3;   % head width
max_length  = 22;
max_width   = 10;

% =============================================
% calculate the arrow head coordinates
% =============================================
den         = x2 - x1 + eps;                                % make sure no devision by zero occurs
teta        = atan( (y2-y1)/den ) + pi*(x2<x1) - pi/2;      % angle of arrow
cs          = cos(teta);                                    % rotation matrix
ss          = sin(teta);
R           = [cs -ss;ss cs];
line_length = sqrt( (y2-y1)^2 + (x2-x1)^2 );                % sizes
head_length = min( line_length*alpha_,max_length );
head_width  = min( line_length*beta,max_length );
x0          = x2*cs + y2*ss;                                % build head coordinats
y0          = -x2*ss + y2*cs;
coords      = R*[x0 x0+head_width/2 x0-head_width/2; y0 y0-head_length y0-head_length];

% =============================================
% plot arrow  (= line + patch of a triangle)
% =============================================
h1          = plot( [x1,x2],[y1,y2], 'Color', [colorin, alphain], 'LineWidth', linewidth);
h2          = plot([coords(1,2), coords(1,1)], [coords(2,2), coords(2,1)], 'Color', [colorin, alphain], 'LineWidth', linewidth);
h3          = plot([coords(1,3), coords(1,1)], [coords(2,3), coords(2,1)], 'Color', [colorin, alphain], 'LineWidth', linewidth);
% h2          = patch( coords(1,:),coords(2,:),[0 0 0] );

end
