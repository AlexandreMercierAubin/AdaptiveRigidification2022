function applyAcceleration(mesh, acc)
    % APPLYACCELERATION applies constant acceleration to the entire mesh
    inds = mesh.ElasticInds;
    if size(mesh.pinnedInds,1) > 0
        triangleRangeIds = 1:mesh.N;
        setDiffIds = false(numel(triangleRangeIds),1);
        setDiffIds(inds) = true;
        setDiffIds = setDiffIds&(1-mesh.pinned);
        inds = triangleRangeIds(setDiffIds);
%         inds = setdiff(inds, mesh.pinned);
    end
    idn = inds * 2;
    ids = idn - 1;

    masses = mesh.mass(idn);
    mesh.f(ids) = mesh.f(ids) + masses * acc(1);
    mesh.f(idn) = mesh.f(idn) + masses * acc(2);

    for i = 1:1:numel(mesh.RigidBodies)
        mesh.RigidBodies(i).Force = mesh.RigidBodies(i).Force + mesh.RigidBodies(i).Mass * acc;
    end
end

