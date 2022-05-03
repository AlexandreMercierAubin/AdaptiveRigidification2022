function [mesh] = loadMeshFromPOLY( fileNameRoot, materials, scale, rotation)
% LOADMESHFROMPOLY Loads POLY ELE and NODE files to create a mesh with the
% given properties
%   fileNameRoot should not include the '.poly' ending.
%   The .ele and .node files MUST exist.
%   Default mesh parameters can be specified or default to the following
%     rho = 100
%     mu, lambda = toLame( 0.48, 2e4 );
%     alpha0, alpha1 = 0, 0
%
%   If loading from poly and constructing the mesh is painfully slow,
%   consider caching the mesh via save( 'filenameroot.mat', mesh ) and load
%   with load( 'filenameroot.mat', mesh );
    if nargin <=1
        materials = [TriangleMaterial()];
    end
    
    if nargin < 3
        scale = [1,1];
    end
    
    if nargin < 4
        rotation = 0;
    end
    
    [ V, I ] = readNODE( sprintf('%s.node', fileNameRoot ) );
    
    V(:,1) = V(:,1) * scale(1);
    V(:,2) = V(:,2) * scale(2);
    
    R = [cosd(rotation), - sind(rotation);...
         sind(rotation), cosd(rotation)];
    
    V = (R*V')';
     
    % I tells us the indices... must be a sequence and start from 1 or we 
    % are in trouble!
    
    [ E, A ] = readELE( sprintf('%s.ele', fileNameRoot ) );

    % NOTE: attributes A are not yet supported, but will be useful for
    % multi-material examples... 

    %[ VV, EE, BME, H ] = readPOLY_triangle( sprintf('%s.poly', fileNameRoot ) ); % H not supported

    % NOTE: poly files could be loaded for the boundary edges, but these
    % are currently found automatically (with some limitations) by the mesh
    % constructor.
    
    mesh = Mesh( V, E, A, materials);

end

