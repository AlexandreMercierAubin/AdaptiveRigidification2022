classdef NullAnimationScript < AnimationScripter
    %NULLANIMATIONSCRIPT use that when there won't be any animation script
    
    properties
        %none
    end
    
    methods
        function obj = NullAnimationScript()
        end
        
        function scriptMesh(obj,mesh3D, integrator, frameNumber, h)
            %default: do nothing
        end
        
        function newDV = scriptVelocity(obj, prevP, prevV, deltaV, activeDOFs, activeDOFsIDs, frame, rigidIDbyVert, h)
            newDV = deltaV;
        end
        
        function newDV = scriptVelocityQuicksolve(obj, prevP, prevV, deltaV, frame, h)
            newDV = deltaV;
        end
        
        function newForce = scriptForceAnimation(obj,force, frame, h)
            newForce = force;
        end
        
        function reset(obj)
        end
    end
end

