function mesh3d = polyLoader(filename, materials, scale, resetMesh, settings)
    %polyLoader(filename, materials, scale, resetMesh, settings)
    %reads .ele and .node files from the .poly format. 
    %If possible reads it from a cache mat file
    if nargin < 4
        resetMesh = false;
    end
    
    if nargin < 3
            scale = [1,1,1];
    end
    scalestring = sprintf('%d_',scale);
    tmpPath2 = ['3d/data/cached/', scalestring, convertStringsToChars(filename),'.mat'];
    if isfile(tmpPath2) && ~resetMesh
        loadedStruct  = load(tmpPath2);
        mesh3d = loadedStruct.mesh3d;
        if isequal(mesh3d.materials, materials)
            return;
        end
    end
    [V,I] = readNODE('3d/data/'+filename+'.node');
    [T,A] = readELE('3d/data/'+filename+'.ele');
    [F,J,K] = boundary_faces(T);
%     tetramesh(T,V);
%     vol = volume(V,T);

    V(:,1) = scale(1)*V(:,1);
    V(:,2) = scale(2)*V(:,2);
    V(:,3) = scale(3)*V(:,3);
    mesh3d = Mesh3D(V,T, A, materials, F, J);
    save(tmpPath2, 'mesh3d');
    settings.recomputeCacheAinv = true;
end