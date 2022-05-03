%assumes the eigenvalues are sorted by eigs
function S = PositiveSemiDefiniteProjection(d2E_div_dsigma2)
    [V, D] = eigs(d2E_div_dsigma2);
    
    if D(1,1) >= 0
        S = d2E_div_dsigma2;
        return;
    end
    
    diagEigenMatrix = D;
    for i = 1:size(diagEigenMatrix,1)
        if diagEigenMatrix(i,i) < 0.0
            diagEigenMatrix(i,i) = 0.0;
        else
            break;
        end
    end
    S = V * diagEigenMatrix / V;
end




