classdef SequentialPositionAnimationScripter < AnimationScripter
    %SequentialPositionAnimationScripter animation scripter that forces
    %nodes to be a specific location at a given frame. The input must be sorted.
    %NOTE that this only works for the last dof of a rigid body
    properties
        frameNumbers = [];
        %assumes sorted DOFs
        dofs = {};% n frame to constrain by constrained node at frame
        positions = {}; %position of node index from nodes at frame
        counter = 1;
        dim = 2;
    end
    
    methods
        function obj = SequentialPositionAnimationScripter()
        end
        
        function obj = scriptMesh(obj, mesh3D, integrator, frame, h)
            isComparison = (obj.counter>=2 && obj.frameNumbers(obj.counter-1) == frame);
            if (obj.counter <= numel(obj.frameNumbers) &&(obj.frameNumbers(obj.counter) == frame) || isComparison)
                frameDOFs = obj.dofs{obj.counter - isComparison};
                mesh3D.animationDOFs = [mesh3D.animationDOFs,frameDOFs];
                if ~isComparison
                    obj.counter = obj.counter + 1;
                end
            end
        end
        
        function newDV = scriptVelocity(obj, prevP, prevV, deltaV, activeDOFs, activeDOFsIDs, frame, rigidIDbyVert, h)
            if frame == 0 
                obj.counter = 1;
            end
            
            newDV = deltaV;
            
            isComparison = (obj.counter>=2 && obj.frameNumbers(obj.counter-1) == frame);
            if (obj.counter <= numel(obj.frameNumbers) &&(obj.frameNumbers(obj.counter) == frame) || isComparison)
                frameDOFs = obj.dofs{obj.counter - isComparison};
                framePositions = obj.positions{obj.counter - isComparison};
%                 DOFsRigid = rigidIDbyVert(frameDOFs);
                velToPos = (framePositions - prevP(frameDOFs))/h;
                
                dvTmp = velToPos - prevV(frameDOFs);
                
                DOFsRigid = rigidIDbyVert(floor((frameDOFs+(obj.dim-1))./obj.dim));
                lastElasticDOF = numel(activeDOFs);
                firstRigid = lastElasticDOF +1;
                
%                 RigidPosChange = zeros(3*max(rigidIDbyVert),1);
                
                for i = 1: numel(DOFsRigid)
                    if DOFsRigid(i) <1 %not rigid
                        newDV(activeDOFs == frameDOFs(i)) = dvTmp(i);
                    else
                        axisDOF = mod(frameDOFs(i)-1,obj.dim)+1;
                        pos = firstRigid + obj.dim * (DOFsRigid(i)-1) + (axisDOF-1);
                        newDV(pos) = dvTmp(i);
                    end
                end
                
            end
        end
        
        function newDV = scriptVelocityQuicksolve(obj, prevP, prevV, deltaV, frame, h)
            newDV = deltaV;
            
            isComparison = (obj.counter>=2 && obj.frameNumbers(obj.counter-1) == frame);
            if (obj.counter <= numel(obj.frameNumbers) &&(obj.frameNumbers(obj.counter) == frame) || isComparison)
                frameDOFs = obj.dofs{obj.counter - isComparison};
                framePositions = obj.positions{obj.counter - isComparison};
                velToPos = (framePositions - prevP(frameDOFs))/h;
                newDV(frameDOFs) = velToPos - prevV(frameDOFs);   
            end
        end
          
        function newForce = scriptForceAnimation(obj,force, frame, h)
            newForce = force;
        end
        
        function reset(obj)
            obj.counter = 1;
        end
    end
end

