function setRigidMotion( mesh, radPerSec, velocity )
%SETRIGIDMOTION Sets mesh point velocities to have given velocity
%   mesh.rigidTransform( radPerSec, velocity )
%   Angular motion is about the center of mass
    totalMass = sum( mesh.mass ) / 2;
    mp = mesh.mass .* mesh.p;
    centerX = sum(mp(1:2:end)) / totalMass;
    centerY = sum(mp(2:2:end)) / totalMass; 
    
    vx = -(mesh.p(2:2:end) - centerY) * radPerSec + velocity(1);
    vy =  (mesh.p(1:2:end) - centerX) * radPerSec + velocity(2);
    mesh.v(1:2:end) = vx;
    mesh.v(2:2:end) = vy;
end