function mesh2d = fetchPoly2D(fileName, resetMesh, materials, scale, rotation, settings)
    %mesh2d = fetchPoly2D(fileName, resetMesh, materials, scale, rotation, settings)
    % If materials have changed, it will recache the data, and change the
    % recomputeCacheAinv in settings to true
    % resetMesh : resets the cache for that mesh
    % materials : list of materials used in the object
    % scale: scaling vector for the scaling in x,y,z
    % rotation: rotation in angles about x,y,z axis
    % settings : the simulation settings, this is passed to toggle
    % recomputedCacheAinv if we detect that it is required. Note that we
    % might not include all the possibilities that require a recaching, so
    % whenever in doubt, delete the cache.
    if nargin < 4
        scale = [1,1];
    end
    
    if nargin < 5
        rotation = 0;
    end
    
    readPath = ['2d/data/' , fileName];
    tmpPath2 = ['2d/data/cached/' , fileName,'_s',string(scale),'_rot',string(rotation),'.mat'];
    tmpPath2 = join(tmpPath2,'');
    if isfile(tmpPath2) && ~resetMesh
        loadedStruct  = load(tmpPath2);
        mesh2d = loadedStruct.mesh2d;
        % check if the mesh matches the defined properties
        if isequal(mesh2d.materials, materials)
            return;
        end
    end
    mesh2d = loadMeshFromPOLY(readPath, materials, scale, rotation);
    save(tmpPath2, 'mesh2d');
    settings.recomputeCacheAinv = true;
end

