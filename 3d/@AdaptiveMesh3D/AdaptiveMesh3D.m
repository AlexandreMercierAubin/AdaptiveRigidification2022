classdef AdaptiveMesh3D < Mesh3D
    properties
        AdaptiveM           % mass matrix with rigid body support
        
        RigidBodies         % array of rigid bodies
        
        ElasticDOFs         % list of elastic DOFs in the 3 * N state vector
        ElasticInds         % indices of elastic vertices

        ElasticTetInds      % indices of elastic elements (duplication ??)
        
        RigidificationValues               % per tet history of E dot values or other metrics, leftmost being most recent
        RigidificationBendValues           % per edge history of dihedral angle
        
        EnableAutomaticRigidification      % setting to enable rigidification
        EnableAutomaticElastification      % setting to enable elastification
        DebugEdotNorms
        VertexDisp                         % local coordinates of rigid dofs in the rigid bodies
        
        isTetElastic                       % list of #triangles booleans representing if the element is elastic or rigidified.
        
        DOFsSpheres
        ActiveDofsCorrespondingID         % list of the activeDofs positions for easier indexing
        ActiveElasticDOFs                 % list of elastic degrees of freedom
        ForceElastifyElementInds          % a list of elements that should not elastify. This is mostly a debugging feature
        
        rigidIDbyVert %output from mexRigidificatior that gives us the rigid body to which a vertex belongs
    end
    
    methods
        function obj = AdaptiveMesh3D(p, t, attributes, materials)
            % ADAPTIVEMESH Extension of Mesh that supports adaptive
            % rigidification. 
            
            if nargin == 1
                % This is as ugly as it gets but I don't see any other way
                % without constructor overloading
                t = 0;
                attributes = [];
                materials = [];
            end
            
            obj@Mesh3D(p, t, attributes, materials);
            
            if isa(p, 'AdaptiveMesh')
                % copy constructor
                fns = properties(p);
                for i = 1:numel(fns)
                    obj.(fns{i}) = p.(fns{i});
                end
                
                for i = 1:numel(obj.RigidBodies)
                    obj.RigidBodies(i) = obj.RigidBodies(i).clone(obj);
                end
                return;
            end
            
            obj.RigidBodies = RigidBody3D.empty;
            obj.ElasticDOFs = 1:obj.N * 3;
            obj.ElasticInds = 1:obj.N;
            obj.ElasticTetInds = 1:size(t, 1);
            obj.AdaptiveM = obj.M;
            
            obj.EnableAutomaticRigidification = 1;
            obj.EnableAutomaticElastification = 1;
            obj.VertexDisp = zeros(obj.N,3);
            obj.isTetElastic = true(size(obj.t,1),1);
            obj.updateRigidState(); % this should not be necessary unless something is not initialized properly but I can't find what is not!
            obj.DOFsSpheres = [];
            obj.ForceElastifyElementInds = [];
            obj.rigidIDbyVert = zeros(obj.N,1);
        end
        
        
        function clone = clone(obj)
            clone = AdaptiveMesh3D(obj);
        end
        
        function M = getM(obj)
            % getM Returns the current mixed elastic rigid mass matrix.
            % Because the rigid bodies have their inertia updated with
            % rotation so that velocities are in the world frame, we need
            % to refresh the rotational inertia blocks in this sparse
            % matrix.
            
            elasticDOFsCount = numel(obj.ElasticDOFs);
            for i = 1:numel(obj.RigidBodies)
                inds = elasticDOFsCount + (i-1)*6 + (4:6);
                obj.AdaptiveM(inds,inds) =  obj.RigidBodies(i).Inertia;
            end
            
            M = obj.AdaptiveM;
        end
        
        function linearMomentum = getAdaptiveLinearMomentum( mesh )
            % GETADAPTIVELINEARMOMENTUM For comparison with the easy
            % computation
            linearMomentum = sum( reshape(mesh.mass(mesh.ElasticDOFs).*mesh.v(mesh.ElasticDOFs), 3, []), 2 );
            for i = 1:numel(mesh.RigidBodies)
                linearMomentum = linearMomentum + mesh.RigidBodies(i).Mass * mesh.RigidBodies(i).Velocity;
            end
        end
        
        function labelElasticTriangles( mesh )
             pt = reshape( mesh.p, 3, mesh.N );
             for i=mesh.ElasticTetInds
                 pos = sum( pt( :, mesh.el(i).t ), 2 ) / 3;
                 text( pos(1), pos(2), pos(3), ""+i );
             end        
        end
        
        function alpha1 = getAlpha1(mesh, elasticOnly)
            alpha1 = [mesh.materials(mesh.materialIndex).alpha1]';
            if elasticOnly
                alpha1 = alpha1(mesh.ElasticTetInds);
                return;
            end
            bigAlpha1 = zeros(size(alpha1,1)*9,1);
            for i = 1:1:size(alpha1,1)
                bigAlpha1(i*9-8:i*9) = alpha1(i);
            end
            alpha1 = bigAlpha1;
        end
        
        function alpha0 = getAlpha0(mesh, elasticOnly)
            if elasticOnly
                alpha0 = mesh.alpha0(mesh.ElasticDOFs);
            else
                alpha0 = mesh.alpha0;
            end
        end
        
        prepare(obj, infMassForPinned);
        
        addForce(mesh, force);
        
        applyAcceleration(obj, acc);
        
        removeRigidBody(obj, body);
        updateParticles(mesh, h, deltav);
        
        [adaptB, gamma] = computeAdaptB(obj);
        
        v = getCurrentVelocity(mesh);
        
        f = getCurrentForce(mesh);
        
        updateRigidState(mesh);
        
        computeRigidForce(mesh);
        
        computeActiveDOFs(mesh);
    end
end

