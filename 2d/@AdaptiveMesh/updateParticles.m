function updateParticles(mesh, h, deltav)
    %UPDATEPARTICLES updates the particles of the mesh given a delta v
    %vector
    elasticN = numel(mesh.ActiveElasticDOFs);
    mesh.v(mesh.ActiveDofsCorrespondingID) = mesh.v(mesh.ActiveDofsCorrespondingID) + deltav(1:elasticN); 
    mesh.p(mesh.ActiveDofsCorrespondingID) = mesh.p(mesh.ActiveDofsCorrespondingID) + h * mesh.v(mesh.ActiveDofsCorrespondingID); 

    index = 1;
    for i = 1:numel(mesh.RigidBodies)
        if mesh.RigidBodies(i).isPinned
            mesh.v(mesh.RigidBodies(i).DOFs) = 0;
            continue;

        end
        mesh.RigidBodies(i).updatePosition(h, deltav(elasticN + index * 3 - 2:elasticN + index * 3));
        index = index + 1;
    end
end