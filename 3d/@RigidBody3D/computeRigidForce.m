function computeRigidForce(body)
    %COMPUTERIGIDFORCE Computes the force and torque on a rigid body using
    %the force of the individual particles
    
    body.Force = body.Force + sum(reshape(body.Mesh.f(body.DOFs), 3, []), 2);
    
    p = reshape(body.Mesh.p(body.DOFs), 3, []);
    torque = p - body.Position;
    f = reshape(body.Mesh.f(body.DOFs), 3, []);
    torque = cross(torque, f);
    
    torqueSum = sum(torque,2);
    body.Torque = body.Torque + torqueSum;
    
    % And also the omega cross J omega term...
    
     
    % Add gyroscopic stabilization term to angular mass
    % http://www8.cs.umu.se/research/reports/2006/005/part1.pdf
    %       L = J omega;
    %       Lhatdt = hat ( L ) * dt;
    %       J += - Lhatdt Jinv Lhat
    % % 	 i.e.,   J += dt*dt*Lhat*Minv*Lhat

    %%  Then.. the torque is    - omega cross Jmod Omega
    %
    % 	massAngular.transform( omega, tmp );
    % 	    	tmp2.cross( tmp, omega );
    % 	    	torque.sub(tmp2);

    % so not the mass we use everywhere else, but the torque we apply here for the - w x J w term.
    
    Jmod = body.Inertia;% - 
    
    body.Torque = body.Torque - cross( body.AngularVelocity, Jmod * body.AngularVelocity );
end
