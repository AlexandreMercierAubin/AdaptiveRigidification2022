function f = getCurrentForce(mesh)
    %get forces on the elatic vertices and rigid bodies
    rigid = zeros(0,0);
    if numel(mesh.RigidBodies) > 0
        rigid = zeros(6*numel(mesh.RigidBodies),1);
        torque = [mesh.RigidBodies.Torque];
        rigid(4:6:end,1) = torque(1,:);
        rigid(5:6:end,1) = torque(2,:);
        rigid(6:6:end,1) = torque(3,:);
        force = [mesh.RigidBodies.Force];
        rigid(1:6:end,1) = force(1,:);
        rigid(2:6:end,1) = force(2,:);
        rigid(3:6:end,1) = force(3,:);
    end
    f = [mesh.f(mesh.ElasticDOFs); rigid];
end

