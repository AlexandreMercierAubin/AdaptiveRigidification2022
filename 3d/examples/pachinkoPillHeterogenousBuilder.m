function baseMesh = pachinkoPillHeterogenousBuilder(materials)
    scale = [1.5,1.5,1.5];
    baseMesh = meshLoader("Pill", [], materials, scale);
    
    [meshSurface,surfaceF,UV,TF,N,NF] = readOBJ('3d/data/Pill.obj');
    meshSurface = meshSurface.*scale;

    [W, Surface2Tet] = Mesh3D.meshSkinning(meshSurface,baseMesh);
    SV = reshape(baseMesh.p',3,[])';
    T = baseMesh.t;
    [F,J,K] = boundary_faces(T);
    baseMesh = Mesh3D(SV,T, [], materials, F, J, W, surfaceF, meshSurface, Surface2Tet);
    baseMesh.setRigidTransform( [0,90,0],[0,0,0]);
    
    v = reshape(baseMesh.p,3,[]);
    n1 = v( :, baseMesh.t(:,1) );
    n2 = v( :, baseMesh.t(:,2) );
    n3 = v( :, baseMesh.t(:,3) );
    n4 = v( :, baseMesh.t(:,4) );
    elCenter = 0.25 * (n1+n2+n3+n4);

    dfromyaxis = elCenter(1,:);
    
    ind = (dfromyaxis > 0);
    attributes = baseMesh.materialIndex;
    attributes(ind) = 2;

    baseMesh.updateMaterials( attributes, materials );
    
    beam0 = Mesh3D(baseMesh);
    beam1 = Mesh3D(baseMesh);
    beam2 = Mesh3D(baseMesh);
    beam3 = Mesh3D(baseMesh);
    beam4 = Mesh3D(baseMesh);
    beam5 = Mesh3D(baseMesh);
    beam6 = Mesh3D(baseMesh);

    minHeight = 23;
    baseMesh.setRigidTransform( [0,0,180],[0.6,-5,minHeight+7]);
    beam0.setRigidTransform( [20,0,0],[0.7,-11,minHeight+4]);
    beam1.setRigidTransform( [0,20,180],[0.4,-13,minHeight+7.5]);
    beam2.setRigidTransform( [0,0,20],[0.5,-9,minHeight+5]);
    beam3.setRigidTransform( [0,0,-20],[0.7,-9,minHeight+3]);
    beam4.setRigidTransform( [-20,-200,0],[0.8,-7,minHeight+2]);
    beam5.setRigidTransform( [90,27,0],[0.8,-3,minHeight+1]);
    beam6.setRigidTransform( [0,0,90],[0.5,-10,minHeight+10]);

    baseMesh.mergeMesh(beam0);
    baseMesh.mergeMesh(beam1);
    baseMesh.mergeMesh(beam2);
    baseMesh.mergeMesh(beam3);
    baseMesh.mergeMesh(beam4);
    baseMesh.mergeMesh(beam5);
    baseMesh.mergeMesh(beam6);
end