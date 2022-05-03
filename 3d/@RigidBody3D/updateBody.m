function updateBody(body, updateMesh)
    %UPDATEBODY updates a rigid body. This must be called after changing
    %the particles a body is composed of. Again, this is expensive so try
    %to minimize the number of calls.

    if nargin < 2
        updateMesh = 1;
    end
    
    if size(body.Mesh.t,2) == 4
        body.TetInds = find(accumarray([body.Mesh.TetsPerParticle{body.Indices}]', 1) == 4)';
    else
        body.TetInds = find(accumarray([body.Mesh.TetsPerParticle{body.Indices}]', 1) == 3)';
    end

    if isempty(body.TetInds)
        if updateMesh
            body.Mesh.removeRigidBody(body);
        else
            body.Mesh.RigidBodies(body.Mesh.RigidBodies == body) = [];
        end
        return;
    end
    
    % this is a fast and good technique to convert any set of indices to dofs
    body.DOFs = reshape([body.Indices * 3 - 2; body.Indices * 3 - 1; body.Indices * 3], 1, []);
    body.isPinned = any(body.Mesh.pinned(body.Indices));
    
    masses = body.Mesh.mass(body.Indices * 3);
    
    body.Mass = sum(masses);
    
    bpx = body.Mesh.p( body.Indices*3-2 );
    bpy = body.Mesh.p( body.Indices*3-1 );
    bpz = body.Mesh.p( body.Indices*3 );
    bvx = body.Mesh.v( body.Indices*3-2 );
    bvy = body.Mesh.v( body.Indices*3-1 );
    bvz = body.Mesh.v( body.Indices*3 );
    
    body.Position = [ dot(bpx, masses); dot(bpy, masses); dot(bpz, masses) ] / body.Mass;
    body.Velocity = [ dot(bvx, masses); dot(bvy, masses); dot(bvz, masses) ] / body.Mass;
    
    % Recall VertexDisp is for the WHOLE mesh... wasteful!! should only
    % keep track of our own vertex displacements... watch out for errors
    % based on this...
    body.Mesh.VertexDisp(body.Indices, 1) = bpx - body.Position(1);
    body.Mesh.VertexDisp(body.Indices, 2) = bpy - body.Position(2);
    body.Mesh.VertexDisp(body.Indices, 3) = bpz - body.Position(3);
    
    bpcom = body.Mesh.VertexDisp( body.Indices, : ); % body particle positions in CoM frame
    
    body.Inertia = zeros(3,3);
    for i = 1: numel(body.Indices)
        rhat = crossProductMatrix( bpcom(i,:) );
        body.Inertia = body.Inertia - masses(i) * rhat^2;
    end
    
    bvcom = zeros( numel(body.Indices), 3 );     % body particle velocities in CoM frame
    bvcom(:,1) = bvx - body.Velocity(1);
    bvcom(:,2) = bvy - body.Velocity(2);
    bvcom(:,3) = bvz - body.Velocity(3);
    
    %rmag = sqrt(sum(bpcom.*bpcom,2));
    
    angularMomentum = zeros( 3, 1 );
    for i = 1: numel(body.Indices)
        angularMomentum = angularMomentum + masses(i) * cross( bpcom(i,:)', bvcom(i,:)' );
    end
    body.AngularVelocity = body.Inertia \ angularMomentum;
    
    body.Inertia0 = body.Inertia;
    body.Rotation = eye(3);
    
    % avoid updating the mesh if not 100% necessary
    if updateMesh
        body.Mesh.updateRigidState();
    end
end 