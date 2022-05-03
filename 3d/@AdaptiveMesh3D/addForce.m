function addForce(mesh, force)
    elasticN2 = numel(mesh.ElasticDOFs);
    mesh.f(mesh.ElasticDOFs) = mesh.f(mesh.ElasticDOFs) + force(1:elasticN2);
    
    for i = 1:numel(mesh.RigidBodies)
        body = mesh.RigidBodies(i);
        body.Force = body.Force + force(elasticN2 + 6 * i - 5:elasticN2 + 6 * i - 3);
        body.Torque = body.Torque + force(elasticN2 + 6 * i - 2: elasticN2 + 6 * i);
    end
end

