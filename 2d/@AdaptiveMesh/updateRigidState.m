function updateRigidState(mesh)    
    % UPDATERIGIDSTATE updates the extended mass matrix to include the
    % current diagonal rigid masses.
    %  Must be called after any alteration to ensure proper
    % behavior of the mesh. Still, this is expensive so avoid updating the
    % mesh for no reason.
   
    %logical setdiff
    nIdRange = [1:mesh.N];
    setDiffIds = true(numel(nIdRange),1);
    setDiffIds( [mesh.RigidBodies.Indices] ) = false;
    mesh.ElasticInds = nIdRange(setDiffIds);
    
    tIdRange = [1:size(mesh.t, 1)];
    setDiffIds = true(numel(tIdRange),1);
    setDiffIds([mesh.RigidBodies.TriInds]) = false;
    mesh.ElasticTriInds = tIdRange(setDiffIds);
    
    mesh.ElasticDOFs = reshape([mesh.ElasticInds * 2 - 1; mesh.ElasticInds * 2], 1, []);
    mesh.ActiveBRows = reshape([mesh.ElasticTriInds * 4 - 3; mesh.ElasticTriInds * 4 - 2; mesh.ElasticTriInds * 4 - 1; mesh.ElasticTriInds * 4], 1, []);
    
    elasticDOFsCount = numel(mesh.ElasticDOFs);
    
    n = elasticDOFsCount + numel(mesh.RigidBodies) * 3;
    
    extendedMdiag = zeros(n, 1);
    extendedMdiag(1:elasticDOFsCount) = mesh.mass(mesh.ElasticDOFs);
    
    if ~isempty(mesh.RigidBodies)
        extendedMdiag(elasticDOFsCount + 1:3:end) = [mesh.RigidBodies.Mass];
        extendedMdiag(elasticDOFsCount + 2:3:end) = [mesh.RigidBodies.Mass];
        extendedMdiag(elasticDOFsCount + 3:3:end) = [mesh.RigidBodies.Inertia];
    end
    
    mesh.AdaptiveM = sparse( 1:n, 1:n, extendedMdiag );
    
    mesh.computeActiveDOFs();
end

