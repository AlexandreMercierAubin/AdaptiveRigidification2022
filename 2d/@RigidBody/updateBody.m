function updateBody(body, updateMesh)
    % UPDATEBODY Updates (recomptues) the mass, inertia, center of mass,
    % position and velocity of a rigid body. This must be called after
    % changing the particles a body is composed of. Again, this is
    % expensive so try to minimize the number of calls.

    if nargin < 2
        updateMesh = 1;
    end
    
    %updateFindAccumarray = tic;
    body.TriInds = find(accumarray([body.Mesh.TrianglesPerParticle{body.Indices}]', 1) == 3)';
    %disp("update accumarray: " + toc(updateFindAccumarray) );
    
    if isempty(body.TriInds)
        if updateMesh
            body.Mesh.removeRigidBody(body);
        else
            body.Mesh.RigidBodies(body.Mesh.RigidBodies == body) = [];
        end
        return;
    end
    
    % this is a fast and good technique to convert any set of indices to dofs
    body.DOFs = reshape([body.Indices * 2 - 1; body.Indices * 2], 1, []);
    body.isPinned = any(body.Mesh.pinned(body.Indices));
    
    masses = body.Mesh.mass(body.Indices * 2);
    
    body.Mass = sum(masses);
    body.Inertia = 0;
    
    %updateposvel = tic;
    body.Position = [dot(body.Mesh.p(body.Indices * 2 - 1), masses); dot(body.Mesh.p(body.Indices * 2), masses)] / body.Mass;
    body.Velocity = [dot(body.Mesh.v(body.Indices * 2 - 1), masses); dot(body.Mesh.v(body.Indices * 2), masses)] / body.Mass;
    %disp("update pos-vel: " + toc(updateposvel) );
    
    % Recall VertexDisp is for the WHOLE mesh... wasteful!! should only
    % keep track of our own vertex displacements... watch out for errors
    % based on this...
    %updatevertex1 = tic;
    body.Mesh.VertexDisp(body.Indices, 1) = body.Mesh.p(body.Indices * 2 - 1) - body.Position(1);
    body.Mesh.VertexDisp(body.Indices, 2) = body.Mesh.p(body.Indices * 2) - body.Position(2);
    %disp("update vertexDisp1: " + toc(updatevertex1) );
    % wastefu but not very costly?
    
    %updatedot = tic;
    body.Inertia = dot(masses, body.Mesh.VertexDisp(body.Indices, 1) .^ 2 + body.Mesh.VertexDisp(body.Indices, 2) .^ 2);
    %disp("update dot: " + toc(updatedot) );  
    % No reason for this to be expensive, no?  Not sure we need to time it
    
    % This looks like the bug!  Should remove the linear velocity part!!
    %  we only want a spin around the center of mass for the angular part!
    %updatevertex2 = tic;
    dxvy = body.Mesh.VertexDisp(body.Indices, 1) .* (body.Mesh.v(body.Indices * 2) - body.Velocity(2) );
    dyvx = body.Mesh.VertexDisp(body.Indices, 2) .* (body.Mesh.v(body.Indices * 2 - 1) - body.Velocity(1) );
    %disp("update vertexDisp2: " + toc(updatevertex2) );
    
    body.Angle = 0;
    body.AngularVelocity = dot(masses, dxvy - dyvx); % cute 2x2 cross product trick
        
    if body.Inertia == 0
        disp('Inerta is zero in updateBody?  This should not happen!');
        body.Inertia = 1;
        body.AngularVelocity = 0;
    else
        body.AngularVelocity = body.AngularVelocity / body.Inertia;
    end
    
    % avoid updating the mesh if not 100% necessary
    %updateRigidStatetime = tic;
    if updateMesh
        body.Mesh.updateRigidState();
    end
    %disp("update rigid state: " + toc(updateRigidStatetime) );
end 