classdef VelocityImpulseAnimationScripter < AnimationScripter
    %VelocityImpulseAnimationScripter animation scripter that create
    %velocities impulses at given frames. The input must be sorted.
    
    properties
        frameNumbers = [];
        dofs = {};% n frame to constrain by constrained node at frame
        velocityImpulse = {}; %position of node index from nodes at frame
        counter = 1;
    end
    
    methods
        function obj = VelocityImpulseAnimationScripter()
        end
        
        function obj = scriptMesh(obj, mesh3D, integrator, frame, h)
            %default: do nothing
            if frame == 0 
                obj.counter = 1;
            end
            isComparison = (obj.counter>=2 && obj.frameNumbers(obj.counter-1) == frame);
            if obj.counter <= numel(obj.frameNumbers) && (obj.frameNumbers(obj.counter) == frame || isComparison)
                frameDOFs = obj.dofs{obj.counter};
                frameImpulse = obj.velocityImpulse{obj.counter};
                mesh3D.v(frameDOFs) = mesh3D.v(frameDOFs)+frameImpulse;
                if ~isComparison
                    obj.counter = obj.counter + 1;
                end
            end
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
            obj.counter = 1;
        end
    end
end

