classdef NeoHookean3DEnergy < EnergyModel
    %CorotationalEnergy class that allows the computation of corotational
    %energy based on https://graphics.pixar.com/library/StableElasticity/paper.pdf
    properties
        
    end
    
    methods
        function obj = NeoHookean3DEnergy()
            obj@EnergyModel();
            obj.name = "neoHookean";
        end 
        
        
        function computeEnergy(obj, mesh3D, deformationGradientF)
            [ii, jj, CblockVals, dpsidF, psi] = mexNeoHookean3D( deformationGradientF, mesh3D.elV, mesh3D.elMu, mesh3D.elLambda );
            C = sparse( ii, jj, CblockVals);
            obj.elasticForces = mesh3D.B' * dpsidF; 
            obj.derivative1Gradient = dpsidF;
            obj.derivative2HessianC = C;
            obj.energy = psi;
        end
        
        function [gradBt,hessBt] = computeIncrementalPotentialGradHess(obj, mu, lambda, x, B)
        %             h2 = h*h;
%             gradBt = h2*obj.elasticForces - 2*fd*h2;
%             hessBt = h2*mesh.B'*(obj.derivative2HessianC*mesh.B);
            Bx = B*x;
            numElements = size(Bx,1)/9;
            
            for i = 1:numElements
                F = reshape(Bx(9*i-8:9*i,:),3,3);

                J = det(F);
                C = F'*F;
                Ic = trace(C);

                %https://graphics.pixar.com/library/StableElasticity/paper.pdf eq 14
                Psi = 0.5 * mu * (Ic-3) + 0.5 * lambda *(J-1)*(J-1);

                gradientPsix = gradient(Psi,x);
                hessianPsix = jacobian(gradientPsix,x);
            end
        end
    end
                      
end

