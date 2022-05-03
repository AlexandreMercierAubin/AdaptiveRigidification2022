function [B, gamma] = getB(obj, cache)
    % COMPUTEADAPTB Computes B for the rigid-elastic combo.
    %   Here, gamma maps to the FULL set of elastic dofs, but we note that
    %   we only actualy use a subset of the dofs to compute B, i.e., those
    %   DOFs that are adjacent to an elastic element, i.e., elements in the
    %   ActiveBRows list.
    
    nelastic = numel(obj.ElasticInds);
    nelastic = nelastic * 3;
    
    nrigid = numel(obj.RigidBodies);
    
    nnz = nelastic + nrigid * 9;
    ii = zeros( nnz, 1 );
    jj = zeros( nnz, 1 );
    vals = zeros( nnz, 1 );

    % This part of the matrix directly copies elastic DOFs ot the full list
    ii(1:nelastic) = obj.ElasticDOFs;
    jj(1:nelastic) = 1:nelastic;
    vals(1:nelastic) = 1;
    nextInd = nelastic+1;
    

    for i = 1:nrigid
        body = obj.RigidBodies(i);
        %[eye(3), cross(-p)] 3 by 6 form matrix
        % ind is y index of elastic element indicies in dofs
        % ind-1 is the x index of elastic elements indices in dofs
        ind3 = body.Indices*3;
        ind2 = ind3-1;
        ind1 = ind3-2;
        % map translational motion of rigid body to its elastic dofs (eye(3))
        count = numel(ind3);
        ii( nextInd:nextInd+count-1 ) = ind1;
        jj( nextInd:nextInd+count-1 ) = nelastic + i*6 - 5;
        vals( nextInd:nextInd+count-1 ) = 1; 
        nextInd = nextInd + count;
        
        ii( nextInd:nextInd+count-1 ) = ind2;
        jj( nextInd:nextInd+count-1 ) = nelastic + i*6 - 4;
        vals( nextInd:nextInd+count-1 ) = 1; 
        nextInd = nextInd + count;
        
        ii( nextInd:nextInd+count-1 ) = ind3;
        jj( nextInd:nextInd+count-1 ) = nelastic + i*6 - 3;
        vals( nextInd:nextInd+count-1 ) = 1; 
        nextInd = nextInd + count;
        
        % map rotational motion of rigid body to its elastic dofs
        %building cross product matrix of -position
       
        p1 = -(body.Mesh.p(ind1) - body.Position(1));
        p2 = -(body.Mesh.p(ind2) - body.Position(2));
        p3 = -(body.Mesh.p(ind3) - body.Position(3));
        %middle top p3
        ii( nextInd:nextInd+count-1 ) = ind2;
        jj( nextInd:nextInd+count-1 ) = nelastic + i*6-2;
        vals( nextInd:nextInd+count-1 ) = p3; 
        nextInd = nextInd + count;
        
        %middle left -p3
        ii( nextInd:nextInd+count-1 ) = ind1;
        jj( nextInd:nextInd+count-1 ) = nelastic + i*6-1;
        vals( nextInd:nextInd+count-1 ) = -p3; 
        nextInd = nextInd + count;
        
        %bottom left p2
        ii( nextInd:nextInd+count-1 ) = ind1;
        jj( nextInd:nextInd+count-1 ) = nelastic + i*6;
        vals( nextInd:nextInd+count-1 ) = p2; 
        nextInd = nextInd + count;
        
        %top right -p2
        ii( nextInd:nextInd+count-1 ) = ind3;
        jj( nextInd:nextInd+count-1 ) = nelastic + i*6-2;
        vals( nextInd:nextInd+count-1 ) = -p2; 
        nextInd = nextInd + count;
        
        %middle right a1
        ii( nextInd:nextInd+count-1 ) = ind3;
        jj( nextInd:nextInd+count-1 ) = nelastic + i*6-1;
        vals( nextInd:nextInd+count-1 ) = p1; 
        nextInd = nextInd + count;
        
        %middle bottom -a1
        ii( nextInd:nextInd+count-1 ) = ind2;
        jj( nextInd:nextInd+count-1 ) = nelastic + i*6;
        vals( nextInd:nextInd+count-1 ) = -p1; 
        nextInd = nextInd + count;
    end
    
    gamma = sparse( ii, jj, vals );
    cache.ActiveB = obj.B(obj.ActiveBRows, :);
    B = cache.ActiveB * gamma;
    
    % I have benchmarked various implementations of the B matrix indexing.
    % In all cases, it ends up being slower than using matlab's builtin
    % indexing. Here are a few examples so no one waste time investigating
    % this further
    
%     dimSelectionm = size(obj.Bt,2);
%     selectionMatrixTransposed = sparse(obj.ActiveBRows,1:numel(obj.ActiveBRows),ones(numel(obj.ActiveBRows),1), dimSelectionm, numel(obj.ActiveBRows));
%     selectionMatrix = sparse(1:numel(obj.ActiveBRows),obj.ActiveBRows,ones(numel(obj.ActiveBRows),1), numel(obj.ActiveBRows), dimSelectionm);
%     gammaTransposed = gamma';
%     if ~isempty(obj.ActiveBRows)
%         B = (gammaTransposed*obj.Bt*selectionMatrixTransposed)';
%         [ii1,jj1,cc1,ii2,jj2,cc2]=mexBSplit(obj.Bt,obj.ActiveBRows');
%     else
%         B=obj.B([], :)*gamma;
%     end
%     tic;B1 = obj.B(obj.ActiveBRows, :) * gamma;toc;
%     tic;B2 = selectionMatrix*obj.B*gamma;toc;
%     tic;B3 = (gammaTransposed*obj.Bt*selectionMatrixTransposed)';toc;
%     tic;B4 = (obj.Bt*selectionMatrixTransposed)'*gamma;toc;
%     tic;B5 = (gammaTransposed*(obj.Bt*selectionMatrixTransposed))';toc;
%     tic;[ii1,jj1,cc1,ii2,jj2,cc2]=mexBSplit(obj.Bt,obj.ActiveBRows'); B6 = sparse(ii1,jj1,cc1,numel(obj.ActiveBRows),size(gamma,1))*gamma;toc;
end

