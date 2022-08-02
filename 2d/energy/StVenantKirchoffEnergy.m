classdef StVenantKirchoffEnergy < EnergyModel
    %StVenantKirchoff3DEnergy class that allows the computation of corotational
    %energy based on Saint-Venant Kirchoff energy formula
    properties
        
    end
    
    methods
        function obj = StVenantKirchoffEnergy()
            obj@EnergyModel();
            obj.name = "StVK2D";
        end 
        
        function computeEnergy(obj, mesh2D, deformationGradientF)
            [ii, jj, CblockVals, dpsidF, psi] = mexComputeSTVKGradHess2D( deformationGradientF, mesh2D.elA, mesh2D.elMu, mesh2D.elLambda );
            sizeC = 4*numel(mesh2D.el);
            C = sparse( ii, jj, CblockVals,sizeC,sizeC);
            obj.elasticForces = mesh2D.B' * dpsidF; 
            obj.derivative1Gradient = dpsidF;
            obj.derivative2HessianC = C;
            obj.energy = psi;
        end
    end
                      
end

