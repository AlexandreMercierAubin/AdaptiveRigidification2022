function baseMesh = forestBuilder(material)
    scale = [1,1.8,1];
    baseMesh = meshLoader("PineTree", [], material, scale, true);
    [meshSurface,surfaceF,UV,TF,N,NF] = readOBJ('3d/data/PineTree.obj');
    meshSurface = meshSurface.*scale;

    [W, Surface2Tet] = Mesh3D.meshSkinning(meshSurface,baseMesh);
    SV = reshape(baseMesh.p',3,[])';
    T = baseMesh.t;
    [F,J,K] = boundary_faces(T);
    baseMesh = Mesh3D(SV,T, [], material(2), F, J, W, surfaceF, meshSurface, Surface2Tet);
    baseMesh.setRigidTransform([90,0,0],[0,0,0]);
    
    %find leaves of pinetree
    v = reshape(baseMesh.p,3,[]);
    n1 = v( :, baseMesh.t(:,1) );
    n2 = v( :, baseMesh.t(:,2) );
    n3 = v( :, baseMesh.t(:,3) );
    n4 = v( :, baseMesh.t(:,4) );
    elCenter = 0.25 * (n1+n2+n3+n4);
    leaves = find(elCenter(3,:)>2.1);
    
    attributes = baseMesh.materialIndex;
    attributes(leaves) = 2; %second in the list of individual mats so tree, leaves

    newMat = [material(2),material(4)];
    baseMesh.updateMaterials( attributes,  newMat);
    
    pine2 = Mesh3D(baseMesh);
    pine3 = Mesh3D(pine2);
    pine4 = Mesh3D(pine2);
    pine8 = Mesh3D(pine2);
    pine9 = Mesh3D(pine2);
   
    scale = [0.5,1.3,0.5];
    baseMesh2 = meshLoader("Tree", [], material, scale, true);
    [meshSurface,surfaceF,UV,TF,N,NF] = readOBJ('3d/data/Tree.obj');
    meshSurface = meshSurface.*scale;

    [W, Surface2Tet] = Mesh3D.meshSkinning(meshSurface,baseMesh2);
    SV = reshape(baseMesh2.p',3,[])';
    T = baseMesh2.t;
    [F,J,K] = boundary_faces(T);
    baseMesh2 = Mesh3D(SV,T, [], material(1), F, J, W, surfaceF, meshSurface, Surface2Tet);
    baseMesh2.setRigidTransform([90,0,0],[0,0,0]);
    
    %find leaves of tree
    v = reshape(baseMesh2.p,3,[]);
    n1 = v( :, baseMesh2.t(:,1) );
    n2 = v( :, baseMesh2.t(:,2) );
    n3 = v( :, baseMesh2.t(:,3) );
    n4 = v( :, baseMesh2.t(:,4) );
    elCenter = 0.25 * (n1+n2+n3+n4);
    leaves = find(elCenter(3,:)>2.5);
    
    attributes = baseMesh2.materialIndex;
    attributes(leaves) = 2;

    newMat = [material(1),material(4)];
    baseMesh2.updateMaterials( attributes,  newMat);
    
    tree2 = Mesh3D(baseMesh2);
    tree3 = Mesh3D(tree2);
    tree5 = Mesh3D(tree2);
    tree6 = Mesh3D(tree2);
    tree7 = Mesh3D(tree2);

    
    pine2.setRigidTransform([0,0,95],[3,-5,0]);
    pine3.setRigidTransform([0,0,175],[-3,-5,0]);
    pine4.setRigidTransform([0,0,29],[0,-10,0]);
    pine8.setRigidTransform([0,0,35],[0,-15,0]);
    pine9.setRigidTransform([0,0,-170],[3,-15,0]);

    baseMesh2.setRigidTransform([0,0,120],[3,0,0.0]);
    tree2.setRigidTransform([0,0,0],[-2.4,0.5,0]);
    tree3.setRigidTransform([0,0,99],[0,-5,0]);
    tree5.setRigidTransform([0,0,-99],[3,-10,0]);
    tree6.setRigidTransform([0,0,129],[-3,-10,0]);
    tree7.setRigidTransform([0,0,60],[-3,-15,0]);

    scale = [0.5,0.5,0.5];
    baseMesh3 = meshLoader("BaseBallBatCoarse", [], material, scale, true);
    [meshSurface,surfaceF,UV,TF,N,NF] = readOBJ('3d/data/BaseballBat.obj');
    meshSurface = meshSurface.*scale;

    [W, Surface2Tet] = Mesh3D.meshSkinning(meshSurface,baseMesh3);
    SV = reshape(baseMesh3.p',3,[])';
    T = baseMesh3.t;
    [F,J,K] = boundary_faces(T);
    baseMesh3 = Mesh3D(SV,T, [], material(3), F, J, W, surfaceF, meshSurface, Surface2Tet);
    baseMesh3.setRigidTransform([0,180,270],[-1.6,2,1.5]);

    baseMesh.mergeMesh(baseMesh2);
    baseMesh.mergeMesh(tree2);
    before = numel(baseMesh.p);
    tbefore = size(baseMesh.t,1);
    baseMesh.mergeMesh(baseMesh3);
    after = numel(baseMesh.p);
    tafter = size(baseMesh.t,1);
    axeRange = before+1:after;
    axeTriRange = tbefore+1:tafter;
    axeRanges = {axeRange,axeTriRange};
    save('3d/data/cached/forestAxeRange','axeRanges');
    baseMesh.mergeMesh(pine2);
    baseMesh.mergeMesh(pine3);
    baseMesh.mergeMesh(tree3);
    baseMesh.mergeMesh(pine4);
    baseMesh.mergeMesh(tree5);
    baseMesh.mergeMesh(tree6);
    baseMesh.mergeMesh(tree7);
    baseMesh.mergeMesh(pine8);
    baseMesh.mergeMesh(pine9);

    zPos = baseMesh.p(3:3:end);
    pinInd = find(zPos(:,1) <= 0.5 );%| zPos(:,1) >= 6.5
    baseMesh.pin(sort(pinInd));
end