function applyAcceleration(mesh, acc)
    % APPLYACCELERATION applies constant acceleration to the entire mesh
    inds = mesh.ElasticInds;
    if size(mesh.pinnedInds,1) > 0
        tetRangeIds = 1:mesh.N;
        setDiffIds = false(numel(tetRangeIds),1);
        setDiffIds(inds) = true;
        setDiffIds = setDiffIds&(1-mesh.pinned);
        inds = tetRangeIds(setDiffIds);
%         inds = setdiff(inds, mesh.pinned);
    end
    id3 = inds * 3;
    id2 = id3 - 1;
    id1 = id3 - 2;

    masses = mesh.mass(id3);
    mesh.f(id1) = mesh.f(id1) + masses * acc(1);
    mesh.f(id2) = mesh.f(id2) + masses * acc(2);
    mesh.f(id3) = mesh.f(id3) + masses * acc(3);
    
    for i = 1:1:numel(mesh.RigidBodies)
        mesh.RigidBodies(i).Force = mesh.RigidBodies(i).Force + mesh.RigidBodies(i).Mass * acc;
    end
end

