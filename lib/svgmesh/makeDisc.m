
r2 = 0.5;

edgeLengthTarget = 0.2;

N2 = ceil( 2*pi*r2 / edgeLengthTarget );

xdata2 = [];
ydata2 = [];
for i=0:(N2-1)
    theta = i/N2*pi*2;
    c = cos(theta)*r2;
    s = sin(theta)*r2;
    xdata2 = [xdata2, c];
    ydata2 = [ydata2, s];
end

shape = polyshape( xdata2, ydata2 );

figure(1);
clf;
subplot(2,1,1);
plot(shape);
hold on;
plot(xdata2,ydata2, '.r' )

axis equal;


filename = "data/disc";

points = [ xdata2; ydata2 ]';
edges2 = [ 1:N2; 2:N2,1 ]; 
edges = [ edges2 ]';
holes = [ ];
% points #V by 2
% edges  #E by 2
% holes  #H by 2
writePOLY_triangle( filename + ".poly", points, edges, holes);

system( pwd + "/lib/svgmesh/triangle.exe  -p -q30 -a0.022 " + filename + ".poly" );
[TV,I] = readNODE( filename + ".1.node" );
[TF,A] = readELE( filename + ".1.ele");

subplot(2,1,2)
if numel(A) > 0
    colours = [ 0.5 0.5 0.5; 1 1 0.5; 0.5 1 1; 1 0.5 1; 1 0.5 0.5; 0.5 0.5 1; 0.5 1 0.5 ];
    faceColor = colours(mod(A,7)+1,:);
    patch('vertices', TV, 'faces', TF, 'edgecol', 'k',  'facecol', 'flat', 'FaceVertexCData', faceColor, 'FaceAlpha', .5, 'EdgeAlpha', .9 );
else
    patch('vertices', TV, 'faces', TF, 'edgecol', 'k',  'facecol', [0.5 0.5 0.5], 'FaceAlpha', .5, 'EdgeAlpha', .9 );

end
axis equal
axis off



