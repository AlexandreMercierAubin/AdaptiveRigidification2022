function [lambda, deltav] = solveLDLTPGS(iterations, Jc, L, D, P, S, b, initial_lambda, cInfo, compliance, td, Winv)
    % SOLVECLDLPGS PGS solver for 2D contact constraints that avoids to
    % build M - h^2K by using a choleski decomposition

    n = size(b, 1);
    deltav = zeros(size(Jc, 2), 1);
    
    lambda = initial_lambda;
    if ~n
        return;
    end 
    

    % What if we build T with rows ONLY for those delta V that affect
    % lambda?  That is, Jc has lots of zero coluns, so we really don't need
    % to know *ALL* dv as we are iterating the soluiton.  In the end, it is
    % as if we really are just building JcAinvJcT... 

%     colWithAllZero = all( Jc == 0 );
%     inds = 1:size(Jc,2);
%     nzix = inds( ~colWithAllZero );
%     JcTnz = Jc(:,nzix )';
    
    % If we are to solve of JcAinvJcT, can we be efficient with the back
    % and forward substitution... we only care about some parts of Ainv
    % so we arrive at a similar problem... 
    
    % perhaps another solution would be to try to approximate JcAinvJcT 
    % with a low rank approximation?  If there are very small values in D, 
    % perhaps it could be natural to truncate these values?
    % Well... we don't have an EVD... both S and D contain scaling and it 
    % is hard to discard values that are not even that small in D without
    % having their importance amplified by S?
    % Likewise, eigenvalue of A don't seem to go down too quickly... :(
    
    % OK.. should try a CG solve for JcAinvJcT
    % is there any opportunity for preconditioning?
        
    
    % Sparse matrix transpose is not free!  probably cheap, but...  
    JcT = Jc';
    
    %JcAinvJcT = Jc * (S * (P * (L' \ (D \ full(L \ (P' * (S * (JcT)))))))) 
    % Or do a CG solve for this?
    
    

    % precompute T = Ainv*Jc' andn ote that this will always be dense!
    %TCompStart = tic;
    T =  S * (P * (L' \ (D \ full(L \ (P' * (S * (JcT)))))));
%     T = Jc\(Jc*T - Winv);
    %td.TComp = toc( TCompStart );

    % Alternatives: an approximate A, where we get a avoid contact
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
%     [ lambda, deltav ] = solveLDLTPGSHelper( iterations, lambda, deltav, T, Dii, b, Jc, mu, compliance );
    [ lambda, deltav ] = mexPGS2D( iterations, lambda, deltav, T, Dii, b, JcT, mu, compliance );
    %td.integrateComputeLambda = toc(ticLambda);

end
