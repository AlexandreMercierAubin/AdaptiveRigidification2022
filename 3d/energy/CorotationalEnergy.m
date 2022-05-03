classdef CorotationalEnergy < EnergyModel
    %CorotationalEnergy class that allows the computation of corotational
    %energy based on %https://graphics.pixar.com/library/StableElasticity/paper.pdf
    properties
        
    end
    
    methods
        function obj = CorotationalEnergy()
            obj@EnergyModel();
            obj.name = "corotational";
        end 
        
        
        function computeEnergy(obj, mesh3D, deformationGradientF)        
            [ii, jj, CblockVals, dpsidF, psi] = mexCorotational3D( deformationGradientF, mesh3D.elV, mesh3D.elMu, mesh3D.elLambda );
            C = sparse( ii, jj, CblockVals);
            obj.elasticForces = mesh3D.B' * dpsidF; 
            obj.derivative1Gradient = dpsidF;
			
			%C = B'*C*B*volume;
            obj.derivative2HessianC = sparse(C);
            obj.energy = psi;
        end
    end
                      
end

