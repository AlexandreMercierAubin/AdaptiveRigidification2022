function removeRigidBody(mesh, body)
    %remove a rigidbody. Perhaps deprecated
    mesh.RigidBodies(mesh.RigidBodies == body) = [];
    mesh.updateRigidState();
end

