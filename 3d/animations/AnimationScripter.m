classdef AnimationScripter < handle
    %ANIMATIONSCRIPTER Class that can force positions/velocities at a specific time
    %step 
    
    properties
        %none needed
    end
    
    methods
        function obj = AnimationScripter()
        end
   end
   methods (Abstract = true)
        scriptMesh(obj,mesh3D, integrator, frameNumber, h);
        newDV = scriptVelocity(obj, prevP, prevV, deltaV, activeDOFs, activeDOFsIDs, frame, rigidIDbyVert, h);
        newDV = scriptVelocityQuicksolve(obj, prevP, prevV, deltaV, frame, h);
        newForce = scriptForceAnimation(obj,force, frame, h);
        reset(obj);
    end
end

