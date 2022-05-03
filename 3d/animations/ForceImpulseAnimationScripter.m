classdef ForceImpulseAnimationScripter < AnimationScripter
    %ForceImpulseAnimationScripter adds force impulses in both the
    %quicksolve and integrator at specific frames
    %frameNumbers, list of frames at which impulses are defined
    %dofs degrees of freedom on which the impulses will act on
    %forceImpulse array of force impulses, should be coupled with dofs (same size and dofs maps to the impulse dofs)
    %counter technical detail that allows me not to iterate through the
    %list of frameNumbers at every frame
    
    properties
        frameNumbers
        dofs
        forceImpulse
        counter
    end
    
    methods
        function obj = ForceImpulseAnimationScripter()
            obj.frameNumbers = [];
            obj.dofs = {};% n frame to constrain by constrained node at frame
            obj.forceImpulse = {}; %position of node index from nodes at frame
            obj.counter = 1;
        end
        
        function scriptMesh(obj,mesh3D, integrator, frameNumber, h)
        end
        
        function newDV = scriptVelocity(obj, prevP, prevV, deltaV, activeDOFs, activeDOFsIDs, frame, rigidIDbyVert, h)
            newDV = deltaV;
        end
        
        function newDV = scriptVelocityQuicksolve(obj, prevP, prevV, deltaV, frame, h)
            newDV = deltaV;
        end
        
        function newForce = scriptForceAnimation(obj,force, frame, h)
            if frame == 0 
                obj.counter = 1;
            end
            
            newForce = force;
            
            isComparison = (obj.counter>=2 && obj.frameNumbers(obj.counter-1) == frame && obj.counter-1 <= numel(obj.frameNumbers));
            if (obj.counter <= numel(obj.frameNumbers) &&(obj.frameNumbers(obj.counter) == frame) || isComparison)
                frameDOFs = obj.dofs{obj.counter - isComparison};
                forceImpulse = reshape(obj.forceImpulse{obj.counter - isComparison},[],1); %doing the reshape just so it does not throw an error is the orientation is not right.
                newForce(frameDOFs') = newForce(frameDOFs') + forceImpulse;
                if ~isComparison
                    obj.counter = obj.counter + 1;
                end
            end
        end
        
        function reset(obj)
            obj.counter = 1;
        end
    end
end

