function computeRigidForce(body)
    %COMPUTERIGIDFORCE Computes the force and torque on a rigid body using
    %the force of the individual particles
    body.Force = body.Force + sum(reshape(body.Mesh.f(body.DOFs), 2, []), 2);

    pxfy = (body.Mesh.p(body.Indices * 2 - 1) - body.Position(1)) .* body.Mesh.f(body.Indices * 2);
    pyfx = (body.Mesh.p(body.Indices * 2)     - body.Position(2)) .* body.Mesh.f(body.Indices * 2 - 1);
    
    body.Torque = sum(pxfy - pyfx); % fast cross product
end
