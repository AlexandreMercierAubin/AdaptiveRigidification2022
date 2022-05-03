function computeActiveDOFs(mesh)

    %COMPUTEACTIVEDOFS Finds dofs that are both unpinned and elastic
    unpinned = mesh.unpinnedDOFs;
    dofIdRange = 1:mesh.N*3;
    setDiffIds = true(numel(dofIdRange),1);
    setDiffIds(mesh.ElasticDOFs) = false;
    rigid = dofIdRange(setDiffIds);
    logical_idx =false(1,numel(mesh.f));
    logical_idx(unpinned) = true;
%     logical_idx(mesh.animationDOFs) = false;
    scaledLogical = logical_idx;
    scaledLogical(rigid) = false;
    logical_idx(rigid) = [];
    
    mesh.ActiveDofsCorrespondingID = scaledLogical;
    mesh.ActiveElasticDOFs = find(logical_idx);
    elasticN = numel(mesh.ElasticDOFs);
    isRigidUnpinned = find(~[mesh.RigidBodies.isPinned]);
    rigidUnpinned = elasticN+reshape([isRigidUnpinned*6-5;isRigidUnpinned*6-4;isRigidUnpinned*6-3;isRigidUnpinned*6-2;isRigidUnpinned*6-1;isRigidUnpinned*6],1,[]);
    mesh.activeDOFs = [mesh.ActiveElasticDOFs, rigidUnpinned];
end