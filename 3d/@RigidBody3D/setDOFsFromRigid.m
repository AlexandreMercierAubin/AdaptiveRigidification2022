function setDOFsFromRigid( body, h )
    %SETDOFSFROMRIGID updates the mesh particles that are part of this body
    %using this body's position and orientation
    
    inds3 = body.Indices * 3;
    inds2 = inds3-1;
    inds1 = inds3-2;
   
    prevP = [body.Mesh.p(inds1),body.Mesh.p(inds2),body.Mesh.p(inds3)]';
    newP = body.Rotation * body.Mesh.VertexDisp(body.Indices, :)' + body.Position;
    body.Mesh.p(inds1) = newP(1,:);
    body.Mesh.p(inds2) = newP(2,:);
    body.Mesh.p(inds3) = newP(3,:);
    
%     omegahat = crossProductMatrix( body.AngularVelocity );

    
    vel = (newP-prevP)/h;
%     vel = omegahat * body.Rotation * body.Mesh.VertexDisp(body.Indices, :)' + body.Velocity;    
    body.Mesh.v(inds1) = vel(1,:);
    body.Mesh.v(inds2) = vel(2,:);
    body.Mesh.v(inds3) = vel(3,:);
end