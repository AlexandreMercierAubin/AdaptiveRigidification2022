classdef RigidBody < handle
    properties
        Mesh            % Mesh of this RigidBody
        
        Position        % rigid position
        Velocity        % rigid body velocity
        Angle           % rigid orientation
        AngularVelocity % rigid angular velocity

        Force           % rigid force accumulator
        Torque          % rigid torque accumulator

        Mass            % mass of rigid
        Inertia         % inertia of the rigid body

        DOFs            % list of rigid DOFs in the 2 * N state vector
        Indices         % indices of rigid vertices 
        TriInds         % rigid triangles indices in this body
        isPinned
    end
    
    methods
        function obj = RigidBody(mesh)
            %RIGIDBODY makes a rigid body to be part of an AdaptiveMesh.
            %   Empty rigid bodies should not exist (rigid bodies with no
            %   triangles)
            obj.Mesh = mesh;
            obj.Position = zeros(2, 1);
            obj.Velocity = zeros(2, 1);
            obj.Angle = zeros(1, 1);
            obj.AngularVelocity = zeros(1, 1);

            obj.Force = zeros(2, 1);
            obj.Torque = zeros(1, 1);

            obj.Mass = 0;
            obj.Inertia = 0; 
            
            obj.DOFs = [];
            obj.Indices = [];
            obj.TriInds = [];  
            obj.isPinned = false;
        end
        
        function clone = clone(obj, mesh)
            clone = RigidBody(mesh);
            fns = properties(obj);
            for i = 1:numel(fns)
                clone.(fns{i}) = obj.(fns{i});
            end
            clone.Mesh = mesh;
        end
        
        function addVertices(obj, indices, updateMesh)
            if nargin < 3
                updateMesh = 1;
            end
            
            obj.Indices = unique([obj.Indices, indices]);
            obj.updateBody(updateMesh);
        end
        
        function removeVertices(obj, indices)
            obj.Indices(ismember(indices, obj.Indices)) = [];
            obj.updateBody();
        end
        
        removeTriangles(obj, tris, updateMesh, settings);
        
        function merge(obj, body, updateMesh)
            if nargin < 3
                updateMesh = 1;
            end
            obj.Mesh.RigidBodies(obj.Mesh.RigidBodies == body) = [];
            obj.Indices = unique([obj.Indices, body.Indices]);
            obj.updateBody(updateMesh);
        end
        
        function updatePosition(body, h, deltav)
            body.Velocity = body.Velocity + deltav(1:2);
            body.AngularVelocity = body.AngularVelocity + deltav(3);            
            body.Position = body.Position + h * body.Velocity;
            body.Angle = body.Angle + h * body.AngularVelocity;
            body.setDOFsFromRigid(h);
        end
        
        function f = computeForceFromRigid(body)
            t = body.Torque;
            r1 = (body.Mesh.p(body.Indices * 2 - 1) - body.Position(1));
            r2 = (body.Mesh.p(body.Indices * 2) - body.Position(2));
            r = (r1.^2 + r2.^2).^0.5;
            force = t./r;
            f = [-force.*(r2./r) + body.Force(1),force.*(r2./r)+ body.Force(2)]./body.Mesh.mass(body.Indices);
        end
        
        %% public function prototypes
        computeRigidForce(body);
        updateBody(obj, updateMesh);
        setDOFsFromRigid(body, h);

    end
end

