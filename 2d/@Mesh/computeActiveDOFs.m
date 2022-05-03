function computeActiveDOFs(mesh)
% active dofs will also have rigid dofs removed (override by adaptive mesh)
    mesh.activeDOFs = mesh.unpinnedDOFs;
end