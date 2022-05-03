function pin(obj, pinnedInds)
    % PIN pins the specified vertex indices
    %   This should perhaps cause an error if a pinned mesh vertex 
    %   is part of a rigid body.
    obj.pinnedInds = reshape([obj.pinnedInds; pinnedInds],[],1);
    obj.pinned = false(obj.N, 1);   % flags pinned indices
    obj.pinned(obj.pinnedInds) = true;
    obj.pinnedDOFs = sort([obj.pinnedInds * 3 - 2; obj.pinnedInds * 3 - 1; obj.pinnedInds * 3]);
    obj.pinnedTets = unique([obj.TetsPerParticle{obj.pinnedInds}]);
    obj.isTetPinned(obj.pinnedTets) = true;
    obj.unpinnedDOFs = setdiff(1:obj.N*3,obj.pinnedDOFs);
    dim = size(obj.t,2);
    pinnedVertPerElement = obj.pinned(obj.t);
    obj.stablePinnedElements = find(sum(pinnedVertPerElement,2) >= dim-1);
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
