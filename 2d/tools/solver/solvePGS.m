function x = solvePGS( iterations, A, b, initial_x, cInfo, h )
    %x = solvePGS( iterations, A, b, initial_x, cInfo )
    %SOLVECPGS PGS solver for 2D contact constraints
    n = size(b, 1);
    
    if nargin > 3
        x = initial_x;
    else
        x = zeros(n, 1);
    end
    
    for it = 1:iterations
        prevX = x;
        for j = 1:n/2
            i = j*2-1;
            delta = dot(A(i, 1:i - 1)', x(1:i - 1)) + dot(A(i, i + 1:n)', x(i + 1:n));
            x(i) = (b(i) - delta) / A(i, i);
            x(i) = max(x(i), 0);
            
            %friction
            k = i+1;
            delta = dot(A(k, 1:k - 1)', x(1:k - 1)) + dot(A(k, k + 1:n)', x(k + 1:n));
            x(k) = (b(k) - delta) / A(k, k);
            hi = cInfo(k/2).frictionCoefficient * x(i);
%             if(-hi >= x(k) || hi <= x(k))
%                 tmp = 0;
%             end
            x(i+1) = min(max(-hi, x(k)), hi);

        end
        if iterations == Inf
            residuals = A*x - b;
            n = numel(b);
            upperBound = zeros(n,1);
            lowerBound = zeros(n,1);
            error = LCP_error(residuals, h*x, diag(A), lowerBound, upperBound);
            if error <= 1e-9 || it > 10000%|| sum(x - prevX) < 1e-30
                break;
            end
        end
    end
end

