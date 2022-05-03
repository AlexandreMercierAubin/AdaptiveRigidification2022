classdef AdaptiveMesh < Mesh
    properties
        AdaptiveM           % mass matrix with rigid body support
        
        RigidBodies         % array of rigid bodies
        
        ElasticDOFs         % list of elastic DOFs in the 2 * N state vector
        ElasticInds         % indices of elastic vertices

        ElasticTriInds      % indices of elastic elements (duplication ??)
        
        RigidificationValues               % per triangle history of E dot values or other metrics, leftmost being most recent
        
        EnableAutomaticRigidification      % setting to enable rigidification
        EnableAutomaticElastification      % setting to enable elastification
        DebugEdotNorms
        VertexDisp      % local coordinates of rigid dofs in the rigid bodies
        isTriElastic    % list of #triangles booleans representing if the element is elastic or rigidified.
        ActiveElasticDOFs   % list of elastic degrees of freedom
        ActiveDofsCorrespondingID %list of the activeDofs positions for easier indexing
        RigidificationDifference %use with difference plotting, shows elements that were either implicitely rigidified or removed from rigidification (causing hinges)
        
        rigidIDbyVert %output from mexRigidificatior that gives us the rigid body to which a vertex belongs
    end
    
    methods
        function obj = AdaptiveMesh(p, t, attributes, materials)
            % ADAPTIVEMESH Extension of Mesh that supports adaptive
            % rigidification. 
            
            if nargin == 1
                % This is as ugly as it gets but I don't see any other way
                % without constructor overloading
                t = 0;
                attributes = [];
                materials = [];
            end
            
            obj@Mesh(p, t, attributes, materials);
            
            if isa(p, 'AdaptiveMesh')
                % copy constructor
                fns = properties(p);
                for i = 1:numel(fns)
                    obj.(fns{i}) = p.(fns{i});
                end
                return;
            end
            obj.VertexDisp = zeros(obj.N,2);
            obj.RigidBodies = RigidBody.empty;
            obj.ElasticDOFs = 1:obj.N * 2;
            obj.ElasticInds = 1:obj.N;
            obj.ActiveElasticDOFs = 1:obj.N*2;
            obj.rigidIDbyVert = zeros(obj.N,1);
            
            obj.ElasticTriInds = 1:size(t, 1);
            obj.AdaptiveM = obj.M;
            
            obj.EnableAutomaticRigidification = 1;
            obj.EnableAutomaticElastification = 1;
            obj.isTriElastic = true(size(obj.t,1),1);
            obj.updateRigidState(); % this should not be necessary unless something is not initialized properly but I can't find what is not!
        end
        
        
        function clone = clone(obj)
            clone = AdaptiveMesh(obj);
        end
        
        function M = getM(obj)
            M = obj.AdaptiveM;
        end
        
        function linearMomentum = getAdaptiveLinearMomentum( mesh )
            % GETADAPTIVELINEARMOMENTUM For comparison with the easy
            % computation
            linearMomentum = sum( reshape(mesh.mass(mesh.ElasticDOFs).*mesh.v(mesh.ElasticDOFs), 2, []), 2 );
            for i = 1:numel(mesh.RigidBodies)
                linearMomentum = linearMomentum + mesh.RigidBodies(i).Mass * mesh.RigidBodies(i).Velocity;
            end
        end
        
        function labelElasticTriangles( mesh )
             pt = reshape( mesh.p, 2, mesh.N );
             for i=mesh.ElasticTriInds
                 pos = sum( pt( :, mesh.el(i).t ), 2 ) / 3;
                 text( pos(1), pos(2), ""+i );
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

