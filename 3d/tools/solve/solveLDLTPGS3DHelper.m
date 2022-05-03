function [lambda, deltav] = solveLDLTPGSHelper( iterations, lambda, deltav, T, Dii, b, Jc, mu, compliance )
    % solveLDLTPGSHelper Performs the inner loop of PGS given friction
    % coefficients.  This version was made to prepare for a mex drop in
    % replacement, but the complexity and number of iterations are low
    % enough that this is not a priority!
    % 
    n = numel(b);
    for it = 1:iterations
        for j = 1:n/3
            i = j*3-2;
            oldLambda = lambda(i);
            lambda(i) = (lambda(i)* Dii(i) - b(i) - Jc(i, :) * deltav) / (Dii(i) + compliance);
            lambda(i) = max(lambda(i), 0);
            deltav = deltav +  T(:,i) * (lambda(i) - oldLambda);
            i = i + 1;
            if mu(j) == 0
                continue
            end
            
            %update lambda for friction
            oldLambda = lambda(i);
            lambda(i) = (lambda(i)* Dii(i) - b(i) - Jc(i, :) * deltav) / (Dii(i) + compliance);
            hi = mu(j) * lambda(i - 1);
            %clamping friction coefficient
            lambda(i) = min(max(-hi, lambda(i)), hi);
            deltav = deltav +  T(:,i) * (lambda(i) - oldLambda);
            
            %update lambda for friction in the second direction
            i = i + 1;
            oldLambda = lambda(i);
            lambda(i) = (lambda(i)* Dii(i) - b(i) - Jc(i, :) * deltav) / (Dii(i) + compliance);
            %clamping friction coefficient
            lambda(i) = min(max(-hi, lambda(i)), hi);
            deltav = deltav +  T(:,i) * (lambda(i) - oldLambda);
        end
    end
end