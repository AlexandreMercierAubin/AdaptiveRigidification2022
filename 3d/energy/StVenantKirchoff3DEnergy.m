classdef StVenantKirchoff3DEnergy < EnergyModel
    %StVenantKirchoff3DEnergy class that allows the computation of corotational
    %energy based on Saint-Venant Kirchoff energy formula
    properties
        
    end
    
    methods
        function obj = StVenantKirchoff3DEnergy()
            obj@EnergyModel();
            obj.name = "StVK";
        end 
        
        function computeEnergy(obj, mesh3D, deformationGradientF)
            [ii, jj, CblockVals, dpsidF, psi] = mexSTVK3D( deformationGradientF, mesh3D.elV, mesh3D.elMu, mesh3D.elLambda );
            C = sparse( ii, jj, CblockVals);
            obj.elasticForces = mesh3D.B' * dpsidF; 
            obj.derivative1Gradient = dpsidF;
            obj.derivative2HessianC = C;
            obj.energy = psi;
        end
    end
                      
end

