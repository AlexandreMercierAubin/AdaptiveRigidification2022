function updateParticles(mesh, h, deltav)
    %UPDATEPARTICLES updates the particles of the mesh given a delta v
    %vector
    elasticN = numel(mesh.ActiveElasticDOFs);

    tmpV = mesh.v(mesh.pinnedDOFs);
    mesh.v(mesh.ActiveDofsCorrespondingID) = mesh.v(mesh.ActiveDofsCorrespondingID) + deltav(1:elasticN); 
    mesh.prevp = mesh.p;
    mesh.p(mesh.ActiveDofsCorrespondingID) = mesh.p(mesh.ActiveDofsCorrespondingID) + h * mesh.v(mesh.ActiveDofsCorrespondingID); 

    index = 1;
    for i = 1:numel(mesh.RigidBodies)
        if mesh.RigidBodies(i).isPinned
            mesh.v(mesh.RigidBodies(i).DOFs) = 0;
            continue;
        end
        mesh.RigidBodies(i).updatePosition(h, deltav(elasticN + index * 6 - 5:elasticN + index * 6));
        index = index + 1;
    end
    if any(tmpV ~= mesh.v(mesh.pinnedDOFs))
        %if this show up, scream!
        print("if this shows up, please scream! Pinned dofs velocity was changed");
        test = tmpV(mesh.pinnedDOFs)
        Noo= mesh.v(mesh.pinnedDOFs)
    end
    tmpP = mesh.p-mesh.prevp;
    if tmpP(mesh.pinnedDOFs) > 0.000001
        print("if this shows up, please scream! Pinned dofs moved");
    end
end