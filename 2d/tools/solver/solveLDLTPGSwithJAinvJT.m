function [lambda, deltav] = solveLDLTPGSwithJAinvJT( iterations, Jc, L, D, P, S, b, initial_lambda, cInfo, compliance )
    % SOLVECLDLPGS PGS solver for 2D contact constraints that uses the LDLT
    % factorization of A = M - h^2K, but does a full assembly of the JcAinvJT matrix 

    n = size(b, 1);
    deltav = zeros(size(Jc, 2), 1);   
    lambda = initial_lambda;
    if ~n
        return;
    end 
    
    JcT = Jc'; % Sparse matrix transpose is not free!  probably cheap, but...  
    JcAinvJcT = Jc * (S * (P * (L' \ (D \ full(L \ (P' * (S * (JcT))))))));
    mu = [ cInfo(:).frictionCoefficient ];
    %[ lambda ] = solveLDLTPGSwithJAinvJTHelper( iterations, JcAinvJcT, lambda, b, mu, compliance );
    [ lambda ] = mexPGS2DwithJAinvJT( iterations, lambda, JcAinvJcT, b, mu, compliance );
    deltav =  S * (P * (L' \ (D \ full(L \ (P' * (S * (JcT * lambda)))))));
end
