classdef EDot3DMexRigidificator < Rigidificator3D
    %EDotMexRigidificator Uses EDot values to tell if a triangle should be
    %made rigid or not. This version uses mex to speedup the process
    properties
        ScaleByMaxEdgeLength
    end
    
    methods
        function obj = EDot3DMexRigidificator()
            %Mexed rigidificator, only uses the elastification call to do
            %both rigidification and elastification in an efficient BFS.
            obj@Rigidificator3D();
        end
        
        function checkForElastification( obj, mesh3D, cache, frame, h, settings )
            if ~isa(mesh3D, 'AdaptiveMesh3D') || ~mesh3D.EnableAutomaticRigidification
                return;
            end
            
            index = mod(frame, obj.FrameCount) + 1;
            %rigidification Edot============================
            % init the Edots
            if ~numel(mesh3D.RigidificationValues)
                    mesh3D.RigidificationValues = zeros(size(mesh3D.t, 1), obj.FrameCount);
            end
            
            if isempty(mesh3D.ElasticTetInds)
                index = mod(frame, obj.FrameCount) + 1;
                mesh3D.RigidificationValues(:, index) = zeros(size(mesh3D.t, 1),1);
                rigidTrisHist = mesh3D.RigidificationValues < obj.RigidificationThreshold;
                cache.edotnorms = [];
            elseif ~settings.RigidificationEnabled
                index = mod(frame, obj.FrameCount) + 1;
                mesh3D.RigidificationValues(:, index) = ones(size(mesh3D.t, 1),1) * obj.RigidificationThreshold;
                rigidTrisHist = mesh3D.RigidificationValues < obj.RigidificationThreshold;
                cache.edotnorms = [];
            else
                index = mod(frame, obj.FrameCount) + 1;
                             
                F = cache.ActiveB * (mesh3D.p + cache.oldp)/2;
                Fdot = cache.ActiveB * mesh3D.v;
                
                EDotNorms = mexEdotNorm3D( F, Fdot );
                cache.edotnorms = EDotNorms; % for plotting
                
                % set EDots
                mesh3D.RigidificationValues(mesh3D.ElasticTetInds, index) = EDotNorms;
                rigidTrisHist = mesh3D.RigidificationValues < obj.RigidificationThreshold;
                if ~settings.RigidificationEnabled
                    rigidTrisHist = mesh3D.RigidificationValues < -inf;
                end
            end
            %End rigidification
            %============================================
            %Elastification Edot ====================
            if isempty(mesh3D.RigidBodies) || ~settings.ElastificationEnabled
                trisToElastify = false(1, size(mesh3D.t,1));
                cache.edotApproxNorms = [];
            else
                rowsIdRange = 1:size(mesh3D.t, 1)* 9;
                setDiffIds = true(numel(rowsIdRange),1);
                setDiffIds(mesh3D.ActiveBRows) = false;
                rigidRows = rowsIdRange(setDiffIds);

                trisIdRange = 1:size(mesh3D.t, 1);
                setDiffIds = setDiffIds(1:9:end);
                rigidTris = trisIdRange(setDiffIds);

                vel = mesh3D.v + cache.ApproximatedDeltaV; 

                FBuffer = mesh3D.B * (mesh3D.p + h*vel+mesh3D.p)/2;
                FApprox = FBuffer(rigidRows);
                FBuffer = mesh3D.B * vel; 
                FdotApprox = FBuffer(rigidRows); 
                
                EDotApproxNorms = mexEdotNorm3D( FApprox, FdotApprox )';
                cache.edotApproxNorms = EDotApproxNorms;

                % rigid tris with value > threshold in Rigidification causes
                % issues, so this prevents that. This might not be necessary
                % though.
                mesh3D.RigidificationValues(rigidTris, index) = 0;
                
                % edot norms again but back into a vector of the entire mesh
                % triangles rather than just the rigids
                globalEDots = zeros(1, size(mesh3D.t, 1));
                globalEDots( rigidTris ) = EDotApproxNorms;

                % new elastic triangles
                trisToElastify = globalEDots > obj.ElastificationThreshold;
                mesh3D.RigidificationValues(trisToElastify, index) = globalEDots(trisToElastify);
            end
            
            %End elastification
            %============================================
                                   
            numRigid = 0;
            if frame >= obj.FrameCount || settings.FirstFrameRigidification
%                 trisToElastify(mesh3D.pinnedTets) = true;%quick test to block pinned tris just in case
                
                if ~isempty(mesh3D.animationInds)
                    animated = false(numel(mesh3D.pinned),1);
                    animated(mesh3D.animationInds) = true;
                    pinnedVertPerTri = mesh3D.pinned(mesh3D.t);
                    stablePinnedElements = find(sum(pinnedVertPerTri,2) >= 3);
                    isStableElement = false(size(mesh3D.t,1),1);
                    isStableElement(stablePinnedElements) = true;
                else
                    stablePinnedElements = mesh3D.stablePinnedElements;
                    isStableElement = false(size(mesh3D.t,1),1);
                    isStableElement(stablePinnedElements) = true;
                    animated = false(numel(mesh3D.pinned),1);
                end
                
                [numRigid, rigidIDbyTet, rigidIDbyVert, isVertElastic, isVertBoundary, isTetElastic, isComponentPinned, isAnimatedComponent] = mexRigidBodyConnectedComponents3D(rigidTrisHist, trisToElastify', mesh3D.AdjagencyMatrix, mesh3D.tint32, mesh3D.N, logical(mesh3D.pinned), mesh3D.valence, stablePinnedElements, isStableElement,animated);
                mesh3D.rigidIDbyVert = rigidIDbyVert;
                assert(sum(isVertElastic) == sum(rigidIDbyVert==0));
                assert(sum(isTetElastic) == sum(rigidIDbyTet <= 0));
                stateChanged = isTetElastic~=mesh3D.isTetElastic;
                mesh3D.RigidificationValues(stateChanged & ~mesh3D.isTetElastic, :) = inf;
                
                hasNotChanged = all(isTetElastic==mesh3D.isTetElastic);
                if hasNotChanged
                    return;
                end
                mesh3D.isTetElastic = isTetElastic;
                
                mesh3D.ElasticTetInds = find(isTetElastic)';
                mesh3D.ElasticInds = find(isVertElastic)';
                rigidVert = find(~isVertElastic);
                
                inds3 = rigidVert * 3;
                RigidDOFs = reshape([inds3 - 2; inds3 - 1; inds3], 1, []);
                
                inds3 = mesh3D.ElasticInds * 3;
                inds9 = mesh3D.ElasticTetInds * 9;
                mesh3D.ElasticDOFs = reshape([inds3 - 2; inds3 - 1; inds3], 1, []);
                mesh3D.ActiveBRows = reshape([inds9 - 8; inds9 - 7; inds9 - 6; inds9 - 5; inds9 - 4; inds9 - 3; inds9 - 2; inds9 - 1; inds9], 1, []);

                [ com, comdot, mass, rotMass, angularMomentum, vertexDisp] = mexRigidBodyProperties3D( numRigid, rigidIDbyVert, mesh3D.p, mesh3D.v, mesh3D.mass );
                
                inds3 = rigidVert * 3;
                mesh3D.VertexDisp(rigidVert,:) = [vertexDisp(inds3-2), vertexDisp(inds3-1), vertexDisp(inds3)];
            end
            
            mesh3D.RigidBodies(numRigid + 1:end) = [];
            lenRigid = numel(mesh3D.RigidBodies);
            for i = 1:numRigid
                if lenRigid < i
                    mesh3D.RigidBodies(i) = RigidBody3D(mesh3D);
                end
                mesh3D.RigidBodies(i).TetInds = find(rigidIDbyTet == i)';
                mesh3D.RigidBodies(i).Indices = find(rigidIDbyVert == i)';
                mesh3D.RigidBodies(i).Position = com(:,i);
                mesh3D.RigidBodies(i).Velocity = comdot(:,i);
                mesh3D.RigidBodies(i).Mass = mass(i);
                inertia = rotMass(:,i);
                mesh3D.RigidBodies(i).Inertia = reshape(inertia,3,3);
                mesh3D.RigidBodies(i).Inertia0 = mesh3D.RigidBodies(i).Inertia;
                mesh3D.RigidBodies(i).AngularVelocity = mesh3D.RigidBodies(i).Inertia \ angularMomentum(:,i);
                mesh3D.RigidBodies(i).Rotation = eye(3);
                
                inds3= mesh3D.RigidBodies(i).Indices * 3;
                mesh3D.RigidBodies(i).DOFs = reshape([inds3 - 2; inds3 - 1; inds3], 1, []);
                mesh3D.RigidBodies(i).Force = [0;0;0];
                mesh3D.RigidBodies(i).Torque = 0;
                mesh3D.RigidBodies(i).isPinned = isComponentPinned(i);
                mesh3D.RigidBodies(i).isAnimated = isAnimatedComponent(i);
            end
           
            
            if frame >= obj.FrameCount || settings.FirstFrameRigidification
                elasticDOFsCount = numel(mesh3D.ElasticDOFs);
                n = elasticDOFsCount + numel(mesh3D.RigidBodies) * 6;
                extendedMdiag = zeros(n, 1);
                extendedMdiag(1:elasticDOFsCount) = mesh3D.mass(mesh3D.ElasticDOFs);

                if ~isempty(mesh3D.RigidBodies)
                    extendedMdiag(elasticDOFsCount + 1:6:end) = [mesh3D.RigidBodies.Mass];
                    extendedMdiag(elasticDOFsCount + 2:6:end) = [mesh3D.RigidBodies.Mass];
                    extendedMdiag(elasticDOFsCount + 3:6:end) = [mesh3D.RigidBodies.Mass];
                end

                mesh3D.AdaptiveM = spdiags(extendedMdiag, 0, n, n);

                for i = 1:numel(mesh3D.RigidBodies)
                    inds = elasticDOFsCount + (i-1)*6 + (4:6);
                    mesh3D.AdaptiveM(inds,inds) =  mesh3D.RigidBodies(i).Inertia;
                end
                
                mesh3D.computeActiveDOFs();
            end
        end
    end
end

