function [R] = rotationFromVectors(v1,v2)
    %ROTATIONFROMVECTORS takes two 3d vectors and generates a rotation matrix
    %that maps from v1 to v2

    nv1 = v1./ norm(v1);
    nv2 = v2./norm(v2);
    crossv = cross(nv1,nv2);
    dotv = dot(nv1,nv2);
    ncrossv = norm(crossv);
    K = [0,-crossv(3), crossv(2);
        crossv(3),0,-crossv(1);
        -crossv(2),crossv(1),0];
    denom = ncrossv*ncrossv;
    if denom > 0
        R = eye(3) + K + K*K*((1-dotv)/(denom));
    else
        R = eye(3);
    end
end

