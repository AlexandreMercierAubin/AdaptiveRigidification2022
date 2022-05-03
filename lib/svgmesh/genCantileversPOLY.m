function [] = genCantileversPOLY()
% GENCANTILEVERSPOLY Utility funciton to run dist mesh to create the points
% and triangles for different cantilever models.
%
%   Examples are saved to the data directory as human readable POLY NODE 
%   and ELE files.

    fd = @(p) drectangle(p, -1.5, 1.5, -1.5, 1.5);
    [p, t] = distmesh2d(fd, @huniform, 0.1, [-1.5, -1.5; 1.5, 1.5], [-1.5, -1.5; 1.5, -1.5; -1.5, 1.5; 1.5, 1.5]);   
    writePT2POLY( '2d/data/square', p, t );
    
    fd = @(p) drectangle(p, -1.5, 3.5, -2, -1.5);
    [p, t] = distmesh2d(fd, @huniform, 0.2, [-1.5, -2; 3.5, -1.5], [-1.5, -2; 3.5, -2; -1.5, -1.5; 3.5, -1.5]);   
    writePT2POLY( '2d/data/barP2', p, t );
    
    fd = @(p) drectangle(p, -1.5, 3.5, -2.5, -1.5);
    [p, t] = distmesh2d(fd, @huniform, 0.2, [-1.5, -2.5; 3.5, -1.5], [-1.5, -2.5; 3.5, -2.5; -1.5, -1.5; 3.5, -1.5]);
    writePT2POLY( '2d/data/longCantileverP2', p, t );
    
    fd = @(p) drectangle(p, -1.5, 1.5, -2.5, -1.5);
    [p, t] = distmesh2d(fd, @huniform, 0.3, [-1.5, -2.5; 1.5, -1.5], [-1.5, -2.5; 1.5, -2.5; -1.5, -1.5; 1.5, -1.5]);
    writePT2POLY('2d/data/cantileverP3', p, t );
       
    fd = @(p) drectangle(p, -1.5, 1.5, -2.5, -1.5);
    [p, t] = distmesh2d(fd, @huniform, 0.2, [-1.5, -2.5; 1.5, -1.5], [-1.5, -2.5; 1.5, -2.5; -1.5, -1.5; 1.5, -1.5]);
    writePT2POLY('2d/data/cantileverP2', p, t);
    
    fd = @(p) drectangle(p, -1.5, 1.5, -2.5, -1.5);
    [p, t] = distmesh2d(fd, @huniform, 0.05, [-1.5, -2.5; 1.5, -1.5], [-1.5, -2.5; 1.5, -2.5; -1.5, -1.5; 1.5, -1.5]);
    writePT2POLY('2d/data/cantileverP05', p, t);
    
    function writePT2POLY( filenameroot, P, T )
        mesh = Mesh( P, T, [],[] ); % create a mesh just for the boundary edges    
        writeNODE( sprintf('%s.node', filenameroot), P );
        writeELE( sprintf('%s.ele', filenameroot), T );  % doesn't support attributes :(
        % POLY file isn't really even used.... but hey!
        edges = vertcat( mesh.boundaryEdges{:} );
        writePOLY_triangle( sprintf('%s.poly', filenameroot), P, edges, [] ); % no holes, and not currently needed, but pehraps later?
    end
end

