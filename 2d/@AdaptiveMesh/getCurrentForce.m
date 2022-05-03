function f = getCurrentForce(mesh)
    %get forces on the elatic vertices and rigid bodies
    rigid = zeros(0,0);
    if size(mesh.RigidBodies,2) > 0
        rigid = zeros(3*size(mesh.RigidBodies,2),1);
        rigid(3:3:end,1) = [mesh.RigidBodies.Torque];
        force = [mesh.RigidBodies.Force];
        rigid(1:3:end,1) = force(1,:);
        rigid(2:3:end,1) = force(2,:);
    end
    f = [mesh.f(mesh.ElasticDOFs); rigid];
end

