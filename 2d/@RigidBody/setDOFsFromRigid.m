function setDOFsFromRigid(body, h)
    %SETDOFSFROMRIGID updates the mesh particles that are part of this body
    %using this body's position and orientation
    rot = [cos(body.Angle), -sin(body.Angle); sin(body.Angle), cos(body.Angle)];

    inds = body.Indices;
    inds2 = inds*2;
    inds1 = inds2 - 1;
    
    pos = body.Position;
    vertexDisp = body.Mesh.VertexDisp(inds, :);
    Rpx = vertexDisp(:, :) * rot(1, :)';
    Rpy = vertexDisp(:, :) * rot(2, :)';

    prevp = [body.Mesh.p(inds1),body.Mesh.p(inds2)]';
    body.Mesh.p(inds1) = pos(1) + Rpx;
    body.Mesh.p(inds2) = pos(2) + Rpy;
    
    vel = ([body.Mesh.p(inds1),body.Mesh.p(inds2)]' - prevp)/h;
    body.Mesh.v(inds1) = vel(1,:);
    body.Mesh.v(inds2) = vel(2,:);
end