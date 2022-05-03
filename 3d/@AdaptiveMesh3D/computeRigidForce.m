function computeRigidForce(mesh)
    for i = 1:1:numel(mesh.RigidBodies)
        mesh.RigidBodies(i).computeRigidForce();
    end
end
