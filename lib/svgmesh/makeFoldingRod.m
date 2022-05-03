function makeFoldingRod()
%MAKEWIGGLEROD 

    offset = 0.2;
    
    amplitude = 0.075;
    frequency = 7;
    halfLength = 2;
    lres = 0.2;
    maxArea = 0.03;    
    filename = "2d/data/foldingCmediumsmall3";

%     offset = 0.04;
%     amplitude = 0.075;
%     frequency = 15;
%     halfLength = 2;
%     lres = 0.02;
%     maxArea = 0.912;    
%     filename = "data/wigglerodlongHD";

    % Not sure what the other parameters were for the short version...
    % or if the geometry is checked in!
    %     lres = 0.039 / 2;
    %     filename = "data/wigglerodHD";

%    f = @(t) [t;amplitude*cos(t*frequency)];
    
    f = @(t) centerline(t);
 
    t = linspace( 0, 3, 200 );
%    t = linspace( -halfLength, halfLength, ceil(halfLength*2/lres) );
    
    c = f(t);

    cx = c(1,:);
    cy = c(2,:);
    
    e = 1e-5;
    % could do symbolic but FD is quick and dirty
    
    t = (f(t+e) - f(t))/e;
    
    tnorm = normrow(t')';
    
    nx = - t(2,:) ./ tnorm;
    ny =   t(1,:) ./ tnorm; 
    
    
    % make an offset curve
    ox1 = cx + offset*nx;
    oy1 = cy + offset*ny;

    ox1 = flip(ox1); % define CCW boundary
    oy1 = flip(oy1);

    ox2 = cx - offset*nx;
    oy2 = cy - offset*ny;
    

    [sox1,soy1] = simplify( ox1, oy1, lres );
    [sox2,soy2] = simplify( ox2, oy2, lres );

    % could also simplify / discretize the 
    % left and right caps that join
    % the upper and lower curve...
    
    capxr = [ sox2(end) sox1(1) ];
    capyr = [ soy2(end) soy1(1) ];

    capxl = [ sox1(end) sox2(1) ];
    capyl = [ soy1(end) soy2(1) ];
    
    [capxr,capyr] = simplify( capxr, capyr, lres );
    [capxl,capyl] = simplify( capxl, capyl, lres );

    % close the curve, avoiding double points at joins
    x = [ sox2(2:end) capxr(2:end) sox1(2:end) capxl(2:end) ];
    y = [ soy2(2:end) capyr(2:end) soy1(2:end) capyl(2:end) ];
    
    ms = 15; 
    figure(1);
    clf
    
    hold off;
    plot( cx, cy, 'b-');
    axis equal
    hold on;

    plot( x,y, 'r.', 'Markersize', ms );
    title( "lres = " + lres );

    % Save this to a poly file, and triangulate!
    N = numel(x);
    points = [ x; y ]';
    edges = [ 1:N; 2:N,1 ]'; 
    holes = [ ];
    % points #V by 2
    % edges  #E by 2
    % holes  #H by 2
    writePOLY_triangle( filename + ".poly", points, edges, holes);    

    system( pwd + "/lib/svgmesh/triangle.exe  -p -q30 -a"+ maxArea + " " + filename + ".poly" );
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
end

function p = centerline( t )
% return the centerline of points t along the curve as a 2 by N matrix
% snake back and forth peridically
%                 \
%                 4
%      /----3-----/
%     2
%      \----1-----
  
%       xl       xh   and r is radius

    xl = -1;
    xh = 1;
    radius = 0.5;

    p = zeros(2,numel(t));
          
    cid = floor(t); % integer curve id
    rid = floor(cid / 4); % repetition id
    sid = mod(cid,4); % segment 0 1 2 or 3
    ct = t - cid;     % curve parameter
    
    ix = sid == 0;
    p(:,ix) =  [ xh*(1-ct(ix)) + xl*ct(ix); rid(ix) * (4*radius) ];

    ix = sid == 1;
    p(:,ix) =  [ -sin(ct(ix)*pi)*radius + xl;     -cos(ct(ix)*pi)*radius + radius + rid(ix) * (4*radius)];
    
    ix = sid == 2;
    p(:,ix) =  [ xl*(1-ct(ix)) + xh*ct(ix); 2*radius + rid(ix) * (4*radius) ];
    
    ix = sid == 3;
    p(:,ix) =  [ sin(ct(ix)*pi)*radius + xh;      -cos(ct(ix)*pi)*radius + 3*radius + rid(ix) * (4*radius)];
    
   
end


function [sx,sy] = simplify( x, y, lres )
% Simplify the given curve to have minimum lres linear resolution but also
% keeping the end points.

% compute total arclength... N = floor of arclengh / lres
% walk along curve making N-1 points between the end points

    el = sqrt( (x(2:end)-x(1:end-1)).^2 + (y(2:end)-y(1:end-1)).^2 ); % edge lengths
    al = sum(el); % arc length
    % max with 1 to make sure we get 2 endpoints!
    N = max( floor( al / lres ), 1 );
    sx = zeros( 1, 2+N-1 );
    sy = zeros( 1, 2+N-1 );
    sx(1) = x(1);
    sy(1) = y(1);
    
    targetlen = al / N;
    
    % There is certainly a a prettier / faster way to do this...
    pal = 0;
    j = 2;
    for i=1:numel(el)
        lx = x(i);
        ly = y(i);
        while pal + el(i) > targetlen
            alpha = (targetlen - pal) / (pal + el(i));
            sx(j) = (1-alpha)*lx+alpha*x(i+1);
            sy(j) = (1-alpha)*ly+alpha*y(i+1);
            lx = sx(j);
            ly = sy(j);
            j = j + 1;
            el(i) = el(i) - (targetlen - pal);                
            pal = 0;
        end
        pal = pal + el(i);
    end

    % clean up the last point
    sx(end) = x(end);
    sy(end) = y(end);

end
            