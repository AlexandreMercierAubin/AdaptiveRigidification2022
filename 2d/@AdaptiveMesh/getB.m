function [B, gamma] = getB(obj, cache)
    % COMPUTEADAPTB Computes B for the rigid-elastic combo.
    %   Here, gamma maps to the FULL set of elastic dofs, but we note that
    %   we only actualy use a subset of the dofs to compute B, i.e., those
    %   DOFs that are adjacent to an elastic element, i.e., elements in the
    %   ActiveBRows list.
    
    nelastic = numel(obj.ElasticDOFs);
    nrigid = numel(obj.RigidBodies);
    
    nnz = nelastic + nrigid * 4;
    ii = zeros( nnz, 1 );
    jj = zeros( nnz, 1 );
    vals = zeros( nnz, 1 );

    % This part of the matrix directly copies elastic DOFs ot the full list
    ii(1:nelastic) = obj.ElasticDOFs;
    jj(1:nelastic) = 1:nelastic;
    vals(1:nelastic) = 1;
    nextInd = nelastic+1;
    
    for i = 1:nrigid
        % ind is y index of elastic element indicies in dofs
        % ind-1 is the x index of elastic elements indices in dofs
        ind = obj.RigidBodies(i).Indices*2;
        % map translational motion of rigid body to its elastic dofs
        count = numel(ind);
        ii( nextInd:nextInd+count-1 ) = ind-1;
        jj( nextInd:nextInd+count-1 ) = nelastic + i*3 - 2;
        vals( nextInd:nextInd+count-1 ) = 1; 
        nextInd = nextInd + count;
        ii( nextInd:nextInd+count-1 ) = ind;
        jj( nextInd:nextInd+count-1 ) = nelastic + i*3 - 1;
        vals( nextInd:nextInd+count-1 ) = 1; 
        nextInd = nextInd + count;
        % map rotational motion of rigid body to its elastic dofs
        ii( nextInd:nextInd+count-1 ) = ind-1;
        jj( nextInd:nextInd+count-1 ) = nelastic + i*3;
        vals( nextInd:nextInd+count-1 ) = - (obj.p(obj.RigidBodies(i).Indices*2) - obj.RigidBodies(i).Position(2)); 
        nextInd = nextInd + count;
        ii( nextInd:nextInd+count-1 ) = ind;
        jj( nextInd:nextInd+count-1 ) = nelastic + i*3;
        vals( nextInd:nextInd+count-1 ) = obj.p(obj.RigidBodies(i).Indices*2-1) - obj.RigidBodies(i).Position(1); 
        nextInd = nextInd + count;        
    end
    
    gamma = sparse( ii, jj, vals );
    cache.ActiveB = obj.B(obj.ActiveBRows, :);
    B = cache.ActiveB * gamma;
end

