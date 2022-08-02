classdef EDiffMexRigidificator < Rigidificator
    % EDiffMexRigidificator Uses E difference norms to tell if an element
    % should be made rigid or not. This version uses mex to speedup the
    % process.
    properties
        Bbuffer
    end
    
    methods
        function obj = EDiffMexRigidificator()
            obj@Rigidificator();
        end
        
        function checkForElastification( obj, mesh2D, cache, frame, h, settings )
            if ~isa(mesh2D, 'AdaptiveMesh') || ~mesh2D.EnableAutomaticRigidification
                return;
            end
            
            index = mod(frame, obj.FrameCount) + 1;
            %rigidification Edot============================
            % init the Edots
            if ~numel(mesh2D.RigidificationValues)
                    mesh2D.RigidificationValues = zeros(size(mesh2D.t, 1), obj.FrameCount);
            end
            
            if isempty(mesh2D.ElasticTriInds)
                index = mod(frame, obj.FrameCount) + 1;
                mesh2D.RigidificationValues(:, index) = zeros(size(mesh2D.t, 1),1);
                rigidTrisHist = mesh2D.RigidificationValues < obj.RigidificationThreshold;
                cache.edotnorms = zeros(size(mesh2D.t,1),1);
            elseif ~settings.RigidificationEnabled
                index = mod(frame, obj.FrameCount) + 1;
                mesh2d.RigidificationValues(:, index) = ones(size(mesh2D.t, 1),1)*obj.RigidificationThreshold;
                rigidTrisHist = mesh2d.RigidificationValues < obj.RigidificationThreshold;
                cache.edotnorms = zeros(size(mesh2D.t,1),1);
            else
                
                Fa = cache.ActiveB * mesh2D.p;
                Fb = cache.ActiveB * cache.oldp; % not so efficient, these would have been previously computed
                
                EDotNorms = mexEdiffNorm2D( Fa, Fb, h );
                
                cache.edotnorms = EDotNorms; % for plotting
                
                % set EDots
                mesh2D.RigidificationValues(mesh2D.ElasticTriInds, index) = EDotNorms;
                rigidTrisHist = mesh2D.RigidificationValues < obj.RigidificationThreshold;
                if ~settings.RigidificationEnabled
                    rigidTrisHist = mesh2D.RigidificationValues < -inf;
                end
            end
            %End rigidification
            %============================================
            %Elastification Edot ====================
            if isempty(mesh2D.RigidBodies) || ~settings.ElastificationEnabled
                trisToElastify = false(1, size(mesh2D.t,1));
            else
                rowsIdRange = 1:size(mesh2D.t, 1) * 4;
                setDiffIds = true(numel(rowsIdRange),1);
                setDiffIds(mesh2D.ActiveBRows) = false;
                rigidRows = rowsIdRange(setDiffIds);

                trisIdRange = 1:size(mesh2D.t, 1);
                setDiffIds = setDiffIds(1:4:end);
                rigidTris = trisIdRange(setDiffIds);

                vel = mesh2D.v + cache.ApproximatedDeltaV; 

                Fa = mesh2D.B * (mesh2D.p + h * vel);
                Fb = mesh2D.B * mesh2D.p; % not so efficient, these would have been previously computed
                Fa = Fa(rigidRows);
                Fb = Fb(rigidRows);
                
                EDotApproxNorms = mexEdiffNorm2D( Fa, Fb, h );

                cache.edotapproxnorms = EDotApproxNorms; % for plotting

                % rigid tris with value > threshold in Rigidification causes
                % issues, so this prevents that. This might not be necessary
                % though.
                mesh2D.RigidificationValues(rigidTris, index) = 0;
                
                % edot norms again but back into a vector of the entire mesh
                % triangles rather than just the rigids
                globalEDots = zeros(1, size(mesh2D.t, 1));
                globalEDots( rigidTris ) = EDotApproxNorms;

                % new elastic triangles
                trisToElastify = globalEDots > obj.ElastificationThreshold;
                
                mesh2D.RigidificationValues(trisToElastify, index) = globalEDots(trisToElastify);
                
            end
            
            %End elastification
            %============================================
                                   
            numRigid = 0;
            if frame >= obj.FrameCount || settings.FirstFrameRigidification
                if obj.PreventPinnedRigidification
                    trisToElastify(mesh2D.pinnedTris) = true;
                end
                if ~isempty(mesh2D.animationInds)                
                    pinned = mesh2D.pinned;
                    pinned(mesh2D.animationInds) = true;
                    pinned = logical(pinned);
                    pinnedVertPerTri = pinned(mesh2D.t);
                    stablePinnedTri = find(sum(pinnedVertPerTri,2) >= 2);
                    isStableTri = false(size(mesh2D.t,1),1);
                    isStableTri(stablePinnedTri) = true;
                else
                    pinned = mesh2D.pinned;
                    pinned = logical(pinned);
                    stablePinnedTri = mesh2D.stablePinnedTri;
                    isStableTri = false(size(mesh2D.t,1),1);
                    isStableTri(stablePinnedTri) = true;
                end
                
                [numRigid, rigidIDbyTri, rigidIDbyVert, isVertElastic, isVertBoundary, isTriElastic, isComponentPinned] = mexRigidConnectedComponents2D(rigidTrisHist, trisToElastify', mesh2D.AdjagencyMatrix, mesh2D.tint32, mesh2D.N , pinned, mesh2D.valence, stablePinnedTri, isStableTri);
                mesh2D.rigidIDbyVert = rigidIDbyVert;
                
                if settings.plotTriImplicitRigidificationElastification
                    preRigidificationDifference = rigidTrisHist(:,1)&rigidTrisHist(:,2)&rigidTrisHist(:,3)&~trisToElastify';
                    mesh2D.RigidificationDifference = preRigidificationDifference - (rigidIDbyTri ~= -1);
                end
                
                assert(sum(isVertElastic) == sum(rigidIDbyVert==0));
                assert(sum(isTriElastic) == sum(rigidIDbyTri<=0));
                assert(all(rigidIDbyTri~=0));
                stateChanged = isTriElastic ~= mesh2D.isTriElastic;
                
                mesh2D.RigidificationValues(stateChanged & ~mesh2D.isTriElastic, :) = inf;
                
                hasNotChanged = all(~stateChanged);
                if hasNotChanged
                    return;
                end
                mesh2D.isTriElastic = isTriElastic;
                
                mesh2D.ElasticTriInds = find(isTriElastic)';
                mesh2D.ElasticInds = find(isVertElastic)';
                
                rigidVert = find(~isVertElastic);
                inds2 = rigidVert * 2;
                RigidDOFs = reshape([inds2 - 1; inds2], 1, []);
                
                rigidVert = find(~isVertElastic);
                
                inds2 = mesh2D.ElasticInds * 2;
                mesh2D.ElasticDOFs = reshape([inds2 - 1; inds2], 1, []);
                
                inds4 = mesh2D.ElasticTriInds * 4;
                mesh2D.ActiveBRows = reshape([inds4 - 3; inds4 - 2; inds4 - 1; inds4], 1, []);
    
                [com, comdot, mass, rotMass, omega, vertexDisp] = mexRigidBodyProperties2D(numRigid, rigidIDbyVert, mesh2D.p, mesh2D.v, mesh2D.mass);
                mesh2D.VertexDisp(rigidVert,2) = vertexDisp(rigidVert*2);
                mesh2D.VertexDisp(rigidVert,1) = vertexDisp(rigidVert*2-1);
            end
            
            mesh2D.RigidBodies(numRigid + 1:end) = [];
            lenRigid = numel(mesh2D.RigidBodies);
            for i = 1:numRigid
                if lenRigid < i
                    mesh2D.RigidBodies(i) = RigidBody(mesh2D);
                end
                mesh2D.RigidBodies(i).TriInds = find(rigidIDbyTri == i)';
                mesh2D.RigidBodies(i).Indices = find(rigidIDbyVert == i)';
                mesh2D.RigidBodies(i).Position = com(:,i);
                mesh2D.RigidBodies(i).Velocity = comdot(:,i);
                mesh2D.RigidBodies(i).Mass = mass(i);
                mesh2D.RigidBodies(i).Inertia = rotMass(i);
                mesh2D.RigidBodies(i).AngularVelocity = omega(i);
                mesh2D.RigidBodies(i).Angle = 0;
                inds2 = mesh2D.RigidBodies(i).Indices * 2;
                mesh2D.RigidBodies(i).DOFs = reshape([inds2 - 1; inds2], 1, []);
                mesh2D.RigidBodies(i).isPinned = isComponentPinned(i);
                mesh2D.RigidBodies(i).Force = [0;0];
                mesh2D.RigidBodies(i).Torque = 0;
            end
           
            if frame >= obj.FrameCount || settings.FirstFrameRigidification
                elasticDOFsCount = numel(mesh2D.ElasticDOFs);
                n = elasticDOFsCount + numel(mesh2D.RigidBodies) * 3;
                extendedMdiag = zeros(n, 1);
                extendedMdiag(1:elasticDOFsCount) = mesh2D.mass(mesh2D.ElasticDOFs);

                if ~isempty(mesh2D.RigidBodies)
                    extendedMdiag(elasticDOFsCount + 1:3:end) = [mesh2D.RigidBodies.Mass];
                    extendedMdiag(elasticDOFsCount + 2:3:end) = [mesh2D.RigidBodies.Mass];
                    extendedMdiag(elasticDOFsCount + 3:3:end) = [mesh2D.RigidBodies.Inertia];
                end

                mesh2D.AdaptiveM = sparse( 1:n, 1:n, extendedMdiag );
                
                mesh2D.computeActiveDOFs();
            end  
        end
    end
end

