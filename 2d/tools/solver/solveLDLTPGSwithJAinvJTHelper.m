function [lambda] = solveLDLTPGSwithJAinvJTHelper( iterations, JAinvJT, lambda, b, mu, complianceIn )
    % solveLDLTPGSHelper Performs the inner loop of PGS given friction
    % coefficients.  This version was made to prepare for a mex drop in
    % replacement.
    % Compliance must be a scalar, or specified per constraint.

    n = numel(b);
    if numel(complianceIn) < n
        compliance = repmat(complianceIn,n,1);
    else
        compliance = complianceIn;
    end
    for it = 1:iterations
        for j = 1:n/2
            i = j*2-1;
            lambda(i) = (lambda(i)* JAinvJT(i,i) - b(i) - JAinvJT(i, :) * lambda) / (JAinvJT(i,i) + compliance(i));
            lambda(i) = max(lambda(i), 0);
            i = i + 1;
            if mu(j) == 0
                continue
            end
            lambda(i) = (lambda(i)* JAinvJT(i,i) - b(i) - JAinvJT(i, :) * lambda) / JAinvJT(i,i);
            hi = mu(j) * lambda(i - 1);
            lambda(i) = min(max(-hi, lambda(i)), hi);
        end
    end
end