classdef NeoHookeanEnergy < EnergyModel
    %CorotationalEnergy class that allows the computation of corotational
    %energy based on https://graphics.pixar.com/library/StableElasticity/paper.pdf
    properties
        
    end
    
    methods
        function obj = NeoHookeanEnergy()
            obj@EnergyModel();
            obj.name = "neoHookean2D";
        end 
        
        
        function computeEnergy(obj, mesh2D, deformationGradientF)
            [ii, jj, CblockVals, dpsidF, psi] = mexNeoHookean2D( deformationGradientF, mesh2D.elA, mesh2D.elMu, mesh2D.elLambda );
            sizeC = 4*numel(mesh2D.el);
            C = sparse( ii, jj, CblockVals,sizeC,sizeC);
            obj.elasticForces = mesh2D.B' * dpsidF; 
            obj.derivative1Gradient = dpsidF;
            obj.derivative2HessianC = C;
            obj.energy = psi;
        end

    end
                      
end

