classdef GravityChangeScript < AnimationScripter
    %GravityChangeScript Stop the simulator from updating the DOF in
    %the gravity direction
    
    properties
        frameNumbers = [];
        gravity = [];
        dofs = {};% n frame to constrain by constrained node at frame
        counter = 1;
    end
    
    methods
        function obj = GravityChangeScript()
        end
        
        function obj = scriptMesh(obj, mesh3D, integrator, frame, h)
            %default: do nothing
            if frame == 0 
                obj.counter = 1;
            end
            isComparison = (obj.counter>=2 && obj.frameNumbers(obj.counter-1) == frame);
            if obj.counter <= numel(obj.frameNumbers) && (obj.frameNumbers(obj.counter) == frame || isComparison)
                integrator.Gravity = obj.gravity(obj.counter);
                obj.counter = obj.counter + 1;
            end
        end
        
        function newDV = scriptVelocity(obj, deltaV, activeDOFsIDs, frameNumber, rigidIDbyVert, h)
            newDV = deltaV;
        end
        
        function newDV = scriptVelocityQuicksolve(obj, deltaV, frameNumber, h)
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

