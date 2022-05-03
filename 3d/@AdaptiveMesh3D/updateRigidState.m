function updateRigidState(mesh)    
    % UPDATERIGIDSTATE updates the extended mass matrix to include the
    % current diagonal rigid masses.
    %  Must be called after any alteration to ensure proper
    % behavior of the mesh. Still, this is expensive so avoid updating the
    % mesh for no reason.
   
    %logical setdiff
    nIdRange = 1:mesh.N;
    setDiffIds = true(numel(nIdRange),1);
    setDiffIds([mesh.RigidBodies.Indices]) = false;
    mesh.ElasticInds = nIdRange(setDiffIds);
    
    tIdRange = 1:size(mesh.t, 1);
    setDiffIds = true(numel(tIdRange),1);
    setDiffIds([mesh.RigidBodies.TetInds]) = false;
    mesh.ElasticTetInds = tIdRange(setDiffIds);
    mesh.isTetElastic = false(size(mesh.t,1),1);
    mesh.isTetElastic(mesh.ElasticTetInds) = true;
    
    dofInds = mesh.ElasticInds * 3; 
    tetInds = mesh.ElasticTetInds * 9;
    mesh.ElasticDOFs = reshape([dofInds - 2; dofInds - 1; dofInds], 1, []);
    mesh.ActiveBRows = reshape([tetInds - 8; tetInds - 7; tetInds - 6; tetInds - 5; tetInds - 4; tetInds - 3; tetInds - 2; tetInds - 1; tetInds], 1, []);
    
    elasticDOFsCount = numel(mesh.ElasticDOFs);
    
    n = elasticDOFsCount + numel(mesh.RigidBodies) * 6;
    
    extendedMdiag = zeros(n, 1);
    extendedMdiag(1:elasticDOFsCount) = mesh.mass(mesh.ElasticDOFs);
    
    if ~isempty(mesh.RigidBodies)
        extendedMdiag(elasticDOFsCount + 1:6:end) = [mesh.RigidBodies.Mass];
        extendedMdiag(elasticDOFsCount + 2:6:end) = [mesh.RigidBodies.Mass];
        extendedMdiag(elasticDOFsCount + 3:6:end) = [mesh.RigidBodies.Mass];
    end
    
    mesh.AdaptiveM = spdiags(extendedMdiag, 0, n, n);
    
    for i = 1:numel(mesh.RigidBodies)
        inds = elasticDOFsCount + (i-1)*6 + (4:6);
        mesh.AdaptiveM(inds,inds) =  mesh.RigidBodies(i).Inertia;
    end
   
    mesh.computeActiveDOFs();
end

