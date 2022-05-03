function mesh3d = fetchMesh3D(fileName, resetMesh, mergeFunction, materials, settings)
    %fetchMesh3D(fileName, resetMesh, mergeFunction, material, settings)
    %returns a mesh, tries to read it from a mat file, otherwise compute it and
    %save it
    md5File = which(func2str(mergeFunction));
    md5Res = mMD5(md5File);
    tmpPath2 = ['3d/data/cached/' , fileName];
    if isfile([tmpPath2,'.mat']) && ~resetMesh
        loadedStruct  = load(tmpPath2);
        mesh3d = loadedStruct.mesh3d;
        if isequal(mesh3d.cacheMergeMaterials, materials) && all(md5Res == mesh3d.cacheMD5)
            return;
        end
    end
    mesh3d = mergeFunction(materials);
    mesh3d.cacheMD5 = md5Res;
    mesh3d.cacheMergeMaterials = materials;
    save(tmpPath2, 'mesh3d');
    settings.recomputeCacheAinv = true;
end

