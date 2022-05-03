function mesh2d = fetchMesh2D(fileName, resetMesh, mergeFunction, materials, settings)
    %fetchMesh2D(fileName, resetMesh, mergeFunction, materials, settings)
    %Looks for a cached mesh in the cached folder. If the mesh is not
    %found, then generate one and save it in the file to save on mesh
    %building the next time we run the scene
    % resetMesh: resets the cache
    % mergeFunction: a builder function to generate a mesh
    % materials: the list of materials used in the scene
    tmpPath2 = ['2d/data/cached/' , fileName];
    
    md5File = which(func2str(mergeFunction));
    md5Res = mMD5(md5File);
    if isfile([tmpPath2,'.mat']) && ~resetMesh
        loadedStruct  = load(tmpPath2);
        mesh2d = loadedStruct.mesh2d;
        
        % check if the mesh materials match the defined properties
        
        if isequal(mesh2d.cacheMergeMaterials, materials) && all(md5Res == mesh2d.cacheMD5)
            return;
        end
    end
    mesh2d = mergeFunction(materials);
    mesh2d.cacheMD5 = md5Res;
    mesh2d.cacheMergeMaterials = materials;
    save(tmpPath2, 'mesh2d');
    settings.recomputeCacheAinv = true;
end

