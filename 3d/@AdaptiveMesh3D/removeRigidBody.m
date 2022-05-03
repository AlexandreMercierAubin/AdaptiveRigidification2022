function removeRigidBody(mesh, body)
% THIS IS NOT CALLED BY ANYONE??
    mesh.RigidBodies(mesh.RigidBodies == body) = [];
    mesh.updateRigidState();
end

