
r1 = 0.75;
r2 = 0.25;

edgeLengthTarget = 0.2;

% the block
xdata = [ -2  2  2 -2 ];
ydata = [ -1 -1  1  1 ];
pos = [ xdata' ydata'];
edges = [ 1:numel(xdata); 2:numel(xdata),1 ]';

holes = 6;
rng(18);
holeix = randperm( 7*3 );
holePosOff = (rand(holes,2) -0.5)*0.4;
holeix = holeix(1:holes);

holeposx =   [ - 4/2 + (ix1(:))/2 + holePosOff(:,1) ];
holeposy =   [ - 2/2 + (ix2(:))/2 + holePosOff(:,2) ];
holepos = [ holeposx, holeposy ];

p = pos;
e = edges;
for j = 1:holes
    r2 = 0.15 + rand(1)*0.15;
    N2 = ceil( 2*pi*r2 / edgeLengthTarget );

    [ix1,ix2] = ind2sub( [7,3], holeix);
    
    cx =  - 4/2 + (ix1(j))/2 + holePosOff(j,1);
    cy =  - 2/2 + (ix2(j))/2 + holePosOff(j,2);
    xdata2 = [];
    ydata2 = [];
    tstart = rand();
    for i=0:(N2-1)
        theta = i/N2*pi*2 + tstart;
        c = cos(theta)*r2 + cx;
        s = sin(theta)*r2 + cy;
        xdata2 = [xdata2, c];
        ydata2 = [ydata2, s];
    end
    pos2 = [xdata2', ydata2'];
    edges2 = [ 1:numel(xdata2); 2:numel(xdata2),1 ]';
    
    [p,e] = combineCurves( p, e, pos2, edges2 );
end

x = [ p(e(:,1),1)'; p(e(:,2),1)'; nan(size(e,1),1)' ];
x = reshape(x,[],1);
y = [ p(e(:,1),2)'; p(e(:,2),2)'; nan(size(e,1),1)' ];
y = reshape(y,[],1);
plot(x,y,'b-');



axis equal;




filename = "2d/data/cheese";


points = p;
edges = e;
holes = holepos;
% points #V by 2
% edges  #E by 2
% holes  #H by 2
writePOLY_triangle( filename + ".poly", points, edges, holes);

system( pwd + "/lib/svgmesh/triangle.exe  -p -q30 -a0.012 " + filename + ".poly" );
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

function [ p, e ] = combineCurves( p1, e1, p2, e2 )
    N = size(p1,1);
    p = [p1;p2];
    e = [e1;e2+N];
end

