classdef CancelDOFAnimationScripter < AnimationScripter
    %CancelDOFAnimationScripter Stop the simulator from updating the DOF
    
    properties
        frameNumbers = [];
        dofs = {};% n frame to constrain by constrained node at frame
        counter = 1;
    end
    
    methods
        function obj = cancelDOFAnimationScript()
        end
        
        function obj = scriptMesh(obj, mesh3D, integrator, frame, h)
            %default: do nothing
            if frame == 0 
                obj.counter = 1;
            end
            isComparison = (obj.counter>=2 && obj.frameNumbers(obj.counter-1) == frame);
            if obj.counter <= numel(obj.frameNumbers) && (obj.frameNumbers(obj.counter) == frame || isComparison)
                frameDOFs = obj.dofs{obj.counter};
                mesh3D.v(frameDOFs) = 0;
                if ~isComparison
                    obj.counter = obj.counter + 1;
                end
                mesh3D.animationDOFs = [mesh3D.animationDOFs,frameDOFs];
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

