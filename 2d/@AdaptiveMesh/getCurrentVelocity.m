function v = getCurrentVelocity(mesh)
    %v = getCurrencyVelocity(mesh)
    % get a stack of elastic velocities and rigid bodies velocities
    rigid = zeros(0,0);
    if numel(mesh.RigidBodies) > 0
        rigid = zeros(3*size(mesh.RigidBodies,2),1);
        rigid(3:3:end,1) = [mesh.RigidBodies.AngularVelocity];
        velocity = [mesh.RigidBodies.Velocity];
        rigid(1:3:end,1) = velocity(1,:);
        rigid(2:3:end,1) = velocity(2,:);
    end
    v = [mesh.v(mesh.ElasticDOFs); rigid];
end

