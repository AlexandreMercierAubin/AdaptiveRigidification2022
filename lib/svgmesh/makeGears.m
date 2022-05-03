
% Note that gear will repeat the first point as the end point to close the
% curve!!!  We will remove that point!!

numTeeth1 = 21;
numTeeth2 = 25;
PA = 15;            % pressure angle (degrees)
Pitch = 2.67;       % teeth pitch (?)
lRes = 0.25;         % linear resolution for discretization
ACDFac = 2;         % Addendum Circle Diameter factor
DCDFac1 = 2.5;      % Dedendum Circle Diameter factor (default = 2)
DCDFac2 = 4;        % Dedendum Circle Diameter factor (default = 2)

[g1x, g1y, PCD1, ACD1] = MakeGear1( numTeeth1, PA, Pitch, lRes, ACDFac, DCDFac1 );
[g2x, g2y, PCD2, ACD2] = MakeGear1( numTeeth2, PA, Pitch, lRes, ACDFac, DCDFac2 );

% NOTE PCD pitch circle diameter is critical for placing gears

% align the gears to mesh with 1 on left and 2 on the right.
% (keep in mind MakeGear1 creates the first tooth centered on top (y).
[g1x, g1y] = rotateZ(g1x, g1y, -pi/2  ); 
[g2x, g2y] = rotateZ(g2x, g2y,  pi/2 + pi/numTeeth2 ); % rotate by 1/2 a tooth to line it up

figure(1);
clf
axis equal;
grid on; hold on;
%plot(g1x, g1y + PCD1/2,'.r');
%plot(g2x, g2y - PCD2/2,'.b');
plot( g1x, g1y ,'-r' );
plot( g2x + PCD2/2 + PCD1/2, g2y, '-b' );

xlim([-15,35]);
ylim([-15,15]);

% Remove the duplicate geometry point before creating geometry
g1x = g1x(1:end-1);
g1y = g1y(1:end-1);
g2x = g2x(1:end-1);
g2y = g2y(1:end-1);


createGearELE( "data/gearA-HD", g1x, g1y, 5, lRes );
createGearELE( "data/gearB-HD", g2x, g2y, 5, lRes );

return

function createGearELE( filename, g1x, g1y, radius, lRes )
    N = length( g1x );
    N1 = ceil( 2*pi*radius / lRes );
    theta = linspace( 2*pi, 0, N1 );
    theta = theta(1:end-1);
    xdata1 = radius*cos(theta);
    ydata1 = radius*sin(theta);
    N1 = numel(theta);

    % combine the edges of the gear and hoe, and make edges
    points = [ g1x, xdata1; g1y, ydata1 ]';
    edges2 = [ 1:N; 2:N,1 ]; 
    edges1 = [ 1:N1; 2:N1,1 ] + [N;N];
    edges = [ edges2, edges1 ]';
    holes = [ 0, 0 ];
    % points #V by 2
    % edges  #E by 2
    % holes  #H by 2
    writePOLY_triangle( filename + ".poly", points, edges, holes);

    system( pwd + "/lib/svgmesh/triangle.exe  -p -q30 -a0.912 " + filename + ".poly" );
    [TV,I] = readNODE( filename + ".1.node" );
    [TF,A] = readELE( filename + ".1.ele");

    figure
    if numel(A) > 0
        colours = [ 0.5 0.5 0.5; 1 1 0.5; 0.5 1 1; 1 0.5 1; 1 0.5 0.5; 0.5 0.5 1; 0.5 1 0.5 ];
        faceColor = colours(mod(A,7)+1,:);
        patch('vertices', TV, 'faces', TF, 'edgecol', 'k',  'facecol', 'flat', 'FaceVertexCData', faceColor, 'FaceAlpha', .5, 'EdgeAlpha', .9 );
    else
        patch('vertices', TV, 'faces', TF, 'edgecol', 'k',  'facecol', [0.5 0.5 0.5], 'FaceAlpha', .5, 'EdgeAlpha', .9 );

    end
    axis equal
    axis off

    % more testing for boundary edges...

%     filename = "data/gearA";
%     [TV,I] = readNODE( filename + ".1.node" );
%     [TF,A] = readELE( filename + ".1.ele");
% 
%     [V,E,BMEo,Ho] = readPOLY_triangleNoVerts( filename + ".poly" );
% 
%     % uh... we should simply take the edges in order that we see them in the
%     boundaries = makeBoundaries(TV,TF);
% 
%     figure(10);
%     clf;
%     for j = 1:numel(boundaries)
%         EE = boundaries{j};
%         for i = 1:size(EE,1);
%             fastQuiver( TV(EE(i,1),:), TV(EE(i,2),:) - TV(EE(i,1),:), [1,0,0], 0.35 ); 
%             hold on;
%         end
%     end
%     axis equal
%     axis off
    
end

