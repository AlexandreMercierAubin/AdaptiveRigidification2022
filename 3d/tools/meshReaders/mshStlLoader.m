function mesh3d = mshStlLoader(filename, attributes, materials,scale, settings)
    %mshStlLoader(filename, attributes, materials, scale)
    %reads a .msh file. If possible reads it from a cache mat file
    if nargin < 5
        resetMesh = false;
    end
    
    if nargin < 4
            scale = [1,1,1];
    end
    scalestring = sprintf('%d_',scale);
    tmpPath2 = ['3d/data/cached/', scalestring, convertStringsToChars(filename),'.mat'];
    if isfile(tmpPath2) && ~resetMesh
        loadedStruct  = load(tmpPath2);
        mesh3d = loadedStruct.mesh3d;
        if isequal(mesh3d.materials, materials) && isequal(mesh3d.materialIndex, attributes)
            return;
        end
    end
    [V,T,F] = readMSH('3d/data/'+filename+'.msh');
    [F,J,K] = boundary_faces(T);

    if nargin < 4
        scale = [1,1,1];
    end
    V(:,1) = scale(1)*V(:,1);
    V(:,2) = scale(2)*V(:,2);
    V(:,3) = scale(3)*V(:,3);
    mesh3d = Mesh3D(V,T, attributes, materials, F, J);
    save(tmpPath2, 'mesh3d');
    settings.recomputeCacheAinv = true;
end