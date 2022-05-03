function v = getCurrentVelocity(mesh)
    rigid = zeros(0,0);
    if numel(mesh.RigidBodies) > 0
        rigid = zeros(6*numel(mesh.RigidBodies),1);
        angular = zeros(3,numel(mesh.RigidBodies));
        for i = 1:numel(mesh.RigidBodies)
            angular(:,i) = [mesh.RigidBodies(i).AngularVelocity];
        end
        rigid(4:6:end,:) = angular(1,:);
        rigid(5:6:end,:) = angular(2,:);
        rigid(6:6:end,:) = angular(3,:);
        velocity = [mesh.RigidBodies.Velocity];
        rigid(1:6:end,:) = velocity(1,:);
        rigid(2:6:end,:) = velocity(2,:);
        rigid(3:6:end,:) = velocity(3,:);
    end
    v = [mesh.v(mesh.ElasticDOFs); rigid];
end

