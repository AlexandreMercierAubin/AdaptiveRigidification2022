function AinvBlocks = computeAInverseDiagonalBlocks(A)
%COMPUTEAINVERSEDIAGONALBLOCKS Computes a sparse matrix with 2x2 blocks of
%only the diagonal of A inverse.
%   There seems to be no really good efficient way to do this.  Probably a
%   good idea to cache the answer.  Could be make slightly faster with a
%   selective forward substitution.

    n = size(A,1);
    
    wb = waitbar(0, 'Computing diagonal blocks of inverse matrix' ); 

    [L, D, P, S] = ldl(A);
    iD = 1 ./ full(diag(D)); % not sure this is faster than backslash??  was this tested?
    
    numBlocks = n/2;    % n must be invisible by 2
    blocks = cell(numBlocks,1);
    for block = 1:numBlocks %parfor this and remove waitbar?
        binds = (block-1)*2+1:block*2;
        b = zeros( n, 2 );
        b(binds,:) = eye(2);
        x = S * (P * (L' \ ( iD .* (L \ (P' * (S * b))))));
        blocks{block} = x(binds,:);
        waitbar( block/numBlocks, wb );
    end
    AinvBlocks = blkdiag( sparse(blocks{1}), blocks{2:end} );
    close(wb);
end

