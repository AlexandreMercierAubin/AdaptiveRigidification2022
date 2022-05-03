function x = solvePGS( iterations, A, b, initial_x, cInfo )
    %solvePGS PGS solver for 3D contact constraints
    n = size(b, 1);
    
    if nargin > 3
        x = initial_x;
    else
        x = zeros(n, 1);
    end
    
    for it = 1:iterations
        for i = 1:n
            sigma = dot(A(i, 1:i - 1)', x(1:i - 1)) + dot(A(i, i + 1:n)', x(i + 1:n));
%             sigma = dot(A(i, :), x) - x(i) * A(i, i);
            x(i) = (b(i) - sigma) / A(i, i);
            
            modi = mod(i, 3);
            if modi == 1 
                x(i) = max(x(i), 0);
            else
                if modi == 2
                    pos = i-1;
                else
                    pos = i-2;
                end
                cinfoPos = floor((i-1)/3) + 1;
                
                hi = cInfo(cinfoPos).frictionCoefficient * x(pos);
                x(i) = min(max(-hi, x(i)), hi);
            end
        end
    end
end

