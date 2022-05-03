function handles = circular_arrow(figHandle, radius, centre, arrow_angle, angle, direction, colour, head_size, head_style)
% This is a function designed to draw a circular arrow onto the current
% figure. It is required that "hold on" must be called before calling this
% function. 
%
% The correct calling syntax is:
%   circular_arrow(height, centre, angle, direction, colour, head_size)
%   where:
%       figHandle - the handle of the figure to be drawn on.
%       radius - the radius of the arrow. 
%       centre - a vector containing the desired centre of the circular
%                   arrow.
%       arrow_angle - the desired orientation angle of the circular arrow.
%                   This is measured in degrees counter-clockwise 
%       angle - the angle between starting and end point of the arrow in
%                   degrees.
%       direction - variable set to determine format of arrow head. Use 1
%                   to get a clockwise arrow, -1 to get a counter clockwise
%                   arrow, 2 to get a double headed arrow and 0 to get just
%                   an arc. 
%       colour (optional) - the desired colour of the arrow, using Matlab's
%                   <a href="matlab:
%                   web('https://au.mathworks.com/help/matlab/ref/colorspec.html')">Color Specification</a>. 
%       head_size (optional) - the size of the arrow head.
%       head_style (optional) - the style of the arrow head.
%                   For more information, see <a href="matlab: 
%                   web('http://au.mathworks.com/help/matlab/ref/annotationarrow-properties.html#property_HeadStyle')">Annotation Arrow Properties</a>.

%Ensure proper number of arguments

%Not sure who made this, it comes from mathworks. Pretty slow to render,
%but nice to use for debugging purpose
if (nargin < 6)||(nargin > 9)
    error(['Wrong number of parameters '...
        'Enter "help circular_arrow" for more information']);
end

% arguments 7, 8 and 9 are optional,
if nargin < 9
   head_style = 'vback2';
end
if nargin < 8
   head_size = 10;
end
if nargin < 7
   colour = 'k';
end

% display a warning if the headstyle has been specified, but direction has
% been set to no heads
if nargin == 9 && direction == 0
    warning(['Head style specified, but direction set to 0! '...
        'This will result in no arrow head being displayed.']);
end

handles = [];

% Check centre is vector with two points
[m,n] = size(centre);
if m*n ~= 2
    error('Centre must be a two element vector');
end


if angle < 0.1
    return;
end

arrow_angle = deg2rad(arrow_angle); % Convert angle to rad
angle = deg2rad(angle); % Convert angle to rad

xc = centre(1);
yc = centre(2);

% Creating (x, y) values that are in the positive direction along the x
% axis and the same height as the centre
x_temp = centre(1) + radius;
y_temp = centre(2);

% Creating x & y values for the start and end points of arc
x1 = (x_temp-xc)*cos(arrow_angle+angle/2) - ...
        (y_temp-yc)*sin(arrow_angle+angle/2) + xc;
x2 = (x_temp-xc)*cos(arrow_angle-angle/2) - ...
        (y_temp-yc)*sin(arrow_angle-angle/2) + xc;
x0 = (x_temp-xc)*cos(arrow_angle) - ...
        (y_temp-yc)*sin(arrow_angle) + xc;
y1 = (x_temp-xc)*sin(arrow_angle+angle/2) + ...
        (y_temp-yc)*cos(arrow_angle+angle/2) + yc;
y2 = (x_temp-xc)*sin(arrow_angle-angle/2) + ... 
        (y_temp-yc)*cos(arrow_angle-angle/2) + yc;
y0 = (x_temp-xc)*sin(arrow_angle) + ... 
        (y_temp-yc)*cos(arrow_angle) + yc;

% Plotting twice to get angles greater than 180
i = 1;

% Creating points
P1 = struct([]);
P2 = struct([]);
P1{1} = [x1;y1]; % Point 1 - 1
P1{2} = [x2;y2]; % Point 1 - 2
P2{1} = [x0;y0]; % Point 2 - 1
P2{2} = [x0;y0]; % Point 2 - 1
centre = [xc;yc]; % guarenteeing centre is the right dimension
n = 1000; % The number of points in the arc
v = struct([]);
    
while i < 3

    v1 = P1{i}-centre;
    v2 = P2{i}-centre;
    c = det([v1,v2]); % "cross product" of v1 and v2
    a = linspace(0,atan2(abs(c),dot(v1,v2)),n); % Angle range
    v3 = [0,-c;c,0]*v1; % v3 lies in plane of v1 and v2 and is orthog. to v1
    v{i} = v1*cos(a)+((norm(v1)/norm(v3))*v3)*sin(a); % Arc, center at (0,0)
    f = plot(v{i}(1,:)+xc,v{i}(2,:)+yc,'Color', colour); % Plot arc, centered at P0
    handles = [handles,f];
    
    i = i + 1;

end

position = struct([]);

% Setting x and y for CW and CCW arrows
if direction == 1
    position{1} = [x2 y2 x2-(v{2}(1,2)+xc) y2-(v{2}(2,2)+yc)];
elseif direction == -1
    position{1} = [x1 y1 x1-(v{1}(1,2)+xc) y1-(v{1}(2,2)+yc)];
elseif direction == 2
    position{1} = [x2 y2 x2-(v{2}(1,2)+xc) y2-(v{2}(2,2)+yc)];
    position{2} = [x1 y1 x1-(v{1}(1,2)+xc) y1-(v{1}(2,2)+yc)];  
elseif direction == 0
    % Do nothing
else
    error('direction flag not 1, -1, 2 or 0.');
end

% Loop for each arrow head
i = 1;
while i < abs(direction) + 1
    h=annotation('arrow'); % arrow head
    set(h,'parent', gca, 'position', position{i}, ...
        'HeadLength', head_size, 'HeadWidth', head_size,...
        'HeadStyle', head_style, 'linestyle','none','Color', colour);
    handles = [handles,h];
    i = i + 1;
end
