function pin(obj, pinnedInds)
    % PIN pins the specified vertex indices
    %   This should perhaps cause an error if a pinned mesh vertex 
    %   is part of a rigid body.
    
    obj.pinnedInds = reshape([obj.pinnedInds; pinnedInds],[],1);
    obj.pinned = false(obj.N, 1);   % flags pinned indices
    obj.pinned(obj.pinnedInds) = true;
    obj.pinnedDOFs = reshape([obj.pinnedInds * 2 - 1; obj.pinnedInds * 2], 1, []);
    obj.pinnedTris = unique([obj.TrianglesPerParticle{obj.pinnedInds}]);
    obj.isTriPinned(obj.pinnedTris) = 1;
    obj.unpinnedDOFs = setdiff(1:obj.N*2, obj.pinnedDOFs);
    pinnedVertPerTri = obj.pinned(obj.t);
    obj.stablePinnedTri = find(sum(pinnedVertPerTri,2) >= 2);
    obj.isStableTri = false(size(obj.t,1),1);
    obj.isStableTri(obj.stablePinnedTri) = true;
    obj.computeActiveDOFs();
    obj.Bii = obj.B(:,obj.unpinnedDOFs);
    obj.Mdii = obj.Md(obj.unpinnedDOFs,obj.unpinnedDOFs);
    obj.Mii = obj.M(obj.unpinnedDOFs,obj.unpinnedDOFs);
    
    for i = 1:numel(obj.objectDOFs)% removes pinnedDOFs from objectDOFsPinned
        dofs = obj.objectDOFsPinned{i};
        pinnedObjectDOFs = ismember(dofs,obj.pinnedDOFs);
        dofs(pinnedObjectDOFs) = [];
        obj.objectDOFsPinned{i} = dofs;
    end
end
