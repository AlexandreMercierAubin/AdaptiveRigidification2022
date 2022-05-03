function precomputeAInverseDiagonalBlocks( cache, mesh, h )
    % precomputeAInverseDiagonalBlocks Builds the rest pose A matrix and 
    % precomputes just the diagonal blocks of the inverse matrix. 
    
    bigB = mesh.B;
    F = mesh.B * mesh.p0;
    
    [ii, jj, CblockVals, ~] = mexComputeSTVKGradHess2D( F, mesh.elA, mesh.elMu, mesh.elLambda );
    bigC = sparse( ii, jj, CblockVals );
    
    a1 = [ mesh.materials(mesh.materialIndex(:)).alpha1 ];
    alpha1 = reshape( repmat( a1, 4, 1 ), [], 1 );
    bigAlpha1 = sparse( 1:numel(alpha1), 1:numel(alpha1), alpha1 );
        
    K  = sparse(bigB' * bigC * bigB);
    Kd = sparse(bigB' * (bigAlpha1 * bigC) * bigB);
    
    M = mesh.M;
    Md = mesh.Md;
       
    A = M - h * (-Md + Kd) - h^2 * K; 

    if isempty(cache.AInvBlocks)
        % this is ONLY needed for boundary nodes... could make faster...
        % Also, at ths point that theis precomp is called, this cache
        % should be empty
        cache.AInvBlocks = computeAInverseDiagonalBlocks( A );
    else
        disp('already have a cached AinvBlocks.. why is precmop getting called?');
    end
    
end