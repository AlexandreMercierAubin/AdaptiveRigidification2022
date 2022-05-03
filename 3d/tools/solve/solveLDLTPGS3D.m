function [lambda, deltav] = solveLDLTPGS3D(iterations, Jc, L, D, P, S, b, initial_lambda, cInfo, compliance, td )
    % SOLVECLDLPGS PGS solver for 2D contact constraints that avoids to
    % build M - h^2K by using a choleski decomposition

    n = size(b, 1);
    deltav = zeros(size(Jc, 2), 1);
    
    lambda = initial_lambda;
    if ~n
        return;
    end 
        
    % Sparse matrix transpose is not free!  probably cheap, but...  
    JcT = Jc';

    % precompute T = Ainv*Jc' and note that this will always be dense!
    %TCompStart = tic;
    T =  S * (P * (L' \ (D \ full(L \ (P' * (S * (JcT)))))));
    %td.TComp = toc( TCompStart );

    % Alternatives: an approximate A, where we get a vaoid contact
    % behaviour but only an approximate dv everywhere else.  Or low rank
    % contact forces... a couple of options.

    %PGSPrepStart = tic;

    % with warm starting, need the deltav to get started!
    deltav =  T * lambda;
    % precompute the values on the diagonal of Jc (M - h^2K)^-1 Jc'
    % Dii must not have any zeros. Make it full for speed and simplicity.
    Dii = full(sum(JcT .* T,1)); 
    
    mu = [ cInfo(:).frictionCoefficient ];
   
    %td.PGSPrep = toc( PGSPrepStart );
        
    %ticLambda = tic;
%     [ lambda, deltav ] = solveLDLTPGS3DHelper( iterations, lambda, deltav, T, Dii, b, Jc, mu, compliance );
    [ lambda, deltav ] = mexPGS3D( iterations, lambda, deltav, T, Dii, b, JcT, mu, compliance );
    %td.integrateComputeLambda = toc(ticLambda);

end
