function addForce(mesh, force)
    % addForce(mesh, force) 
    % function of mesh that adds a force array to the particles including
    % rigid bodies
    % the force array should be of size 2*elasticElements+3*rigidBodies
    elasticN2 = numel(mesh.ElasticDOFs);
    mesh.f(mesh.ElasticDOFs) = mesh.f(mesh.ElasticDOFs) + force(1:elasticN2);
    
    for i = 1:numel(mesh.RigidBodies)
        body = mesh.RigidBodies(i);
        body.Force = body.Force + force(elasticN2 + 3 * i - 2:elasticN2 + 3 * i - 1);
        body.Torque = body.Torque + force(elasticN2 + 3 * i);
    end
end

