function baseMesh = octopusBowlBuilder(material)
    scale = [10,10,10];
    baseMesh = meshLoader("octopusP2", [], material, scale);
    [meshSurface,surfaceF,UV,TF,N,NF] = readOBJ('3d/data/octopus.obj');
    meshSurface = meshSurface.*scale;

    [W, Surface2Tet] = Mesh3D.meshSkinning(meshSurface,baseMesh);
    SV = reshape(baseMesh.p',3,[])';
    T = baseMesh.t;
    [F,J,K] = boundary_faces(T);
    baseMesh = Mesh3D(SV,T, [], material(1), F, J, W, surfaceF, meshSurface, Surface2Tet);

    baseMesh.setRigidTransform([90,0,0],[0,-1,5],true);
end