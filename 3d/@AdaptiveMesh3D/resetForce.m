function resetForce(mesh)
    %reset the forces to zero for elastic dofs and rigid bodies
    mesh.f(:) = 0;
    
    for i = 1:numel(mesh.RigidBodies)
        mesh.RigidBodies(i).Force = [0; 0; 0];
        mesh.RigidBodies(i).Torque = [0; 0; 0];
    end
end

