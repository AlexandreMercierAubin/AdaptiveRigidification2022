function x = PCGnoLoop( A, b, cache, blockDOFs )
% PCGNOLOOP One iteration of Conjugate Gradient on Ax=b, using a
% preconditioner and starting with the initial point (should be zero! Use
% shifted Krylov to warm start instead).
%   A           matlab function that multiplies by A
%   b           the right hand side
%   cache       has the preconditioner, a matlab function that applies the 
%               preconditioner. T
    x = zeros(size(b));    
    r = b; % first x is zero, so residual is b
    z = cache.preconditioner(r);
    Az = A(z);
    
    prevN = 0;
    for i = 1:numel(blockDOFs)
        %here we are using numel in case some DOFs of the full system are
        %pinned
        N = numel(blockDOFs{i});
        dofs = (prevN+1:prevN+N)';
        zAz = (z(dofs)'*Az(dofs));
        if zAz == 0
            x(dofs) = x(dofs) + 0;
        else
            x(dofs) = x(dofs) + ((r(dofs)' * z(dofs)) / zAz) * z(dofs);
        end
        prevN = prevN+N;
    end
end