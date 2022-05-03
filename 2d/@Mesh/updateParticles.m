function updateParticles(mesh, h, deltav)
    %updateParticles(mesh, h, deltav) advances positions of the particles in time using the change in
    %velocities provided
    mesh.v(mesh.activeDOFs) = mesh.v(mesh.activeDOFs) + deltav;
    mesh.v(mesh.pinnedDOFs) = mesh.v(mesh.pinnedDOFs) * 0;
    mesh.p = mesh.p + h * mesh.v;
end