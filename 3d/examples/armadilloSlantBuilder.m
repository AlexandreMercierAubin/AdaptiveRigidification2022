function baseMesh = armadilloSlantBuilder(material)

    armadilloFrame = polyLoader("arma_6", material(1),[1.0,1.00,1.0]);
    [armadilloSurface,surfaceF,UV,TF,N,NF] = readOBJ('3d/data/armadillo.obj');

    [W, Surface2Tet] = Mesh3D.meshSkinning(armadilloSurface,armadilloFrame);
    SV = reshape(armadilloFrame.p',3,[])';
    T = armadilloFrame.t;
    [F,J,K] = boundary_faces(T);
    baseMeshTmp = Mesh3D(SV,T, [], material(1), F, J, W, surfaceF, armadilloSurface, Surface2Tet);

    baseMesh2 = Mesh3D(baseMeshTmp);
    baseMeshTmp.setRigidTransform([270,0,90],[0,0,1.5]);
    baseMesh2.setRigidTransform([10,0,30],[0,0,1.5]);

    baseMesh = baseMeshTmp;
    baseMesh.mergeMesh(baseMesh2);
end