function setRigidMotion( mesh, omega, velocity )
%SETRIGIDMOTION Sets mesh point velocities to have given velocity
%   mesh.rigidTransform( omega, velocity )
%   Angular motion is about the center of mass
%   Both omega and v and 3 component vectors in the global frame

    if ( size(omega,1) == 1 )
        omega = omega';
    end
    if ( size(velocity,1) ==1 )
        velocity = velocity';
    end

    totalMass = sum( mesh.mass ) / 3;
    mp = mesh.mass .* mesh.p;
    centerPos = zeros(3,1);
    centerPos(1) = sum(mp(1:3:end)) / totalMass;
    centerPos(2) = sum(mp(2:3:end)) / totalMass; 
    centerPos(3) = sum(mp(3:3:end)) / totalMass; 

    for i = 1:mesh.N
        ind = i*3-2:i*3;
        mesh.v(ind) = cross( omega, mesh.p(ind) - centerPos ) + velocity;
    end    
end