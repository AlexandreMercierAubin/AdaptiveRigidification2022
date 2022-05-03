function [lambda, deltav] = solveLDLTPGSHelper( iterations, lambda, deltav, T, Dii, b, Jc, mu, complianceIn )
    % solveLDLTPGSHelper Performs the inner loop of PGS given friction
    % coefficients.  This version was made to prepare for a mex drop in
    % replacement, but the complexity and number of iterations are low
    % enough that this is not a priority!
    initialLambda = lambda;
    n = numel(b);
    if numel(complianceIn) < n
        compliance = repmat(complianceIn,n,1);
        compliance = diag(compliance);
    else
        compliance = complianceIn;
    end
    for it = 1:iterations
        for j = 1:n/2
            i = j*2-1;
            oldLambda = lambda(i);
            lambda(i) = (lambda(i)* Dii(i) - b(i) - Jc(i, :) * deltav) / (Dii(i) + compliance(i,i));
            lambda(i) = max(lambda(i), 0);
            deltav = deltav +  T(:,i) * (lambda(i) - oldLambda);
            i = i + 1;
            if mu(j) == 0
                continue
            end
            oldLambda = lambda(i);
            lambda(i) = (lambda(i)* Dii(i) - b(i) - Jc(i, :) * deltav) / (Dii(i) + compliance(i,i));
            hi = mu(j) * lambda(i - 1);
            lambda(i) = min(max(-hi, lambda(i)), hi);
            deltav = deltav +  T(:,i) * (lambda(i) - oldLambda);
        end
    end
end