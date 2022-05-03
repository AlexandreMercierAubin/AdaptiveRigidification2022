function [angle] = dihedralAngle(v1,v2,o1,o2)
%DIHEDRALANGLE compute the dihedralAngle from two triangles separated by an
%edge. v1 & v2 vertices of edge. o1 & o2 vertices opposite of edge 
    b1 = o1 - v1;
    b2 = v2 - v1;
    b3 = o2 - v1;
    n1 = cross(b1,b2);
    n1 = n1./norm(n1);
    
    n2 = cross(b2,b3);
    n2 = n2./norm(n2);
    
    e = b2./norm(b2);
    m1 = cross(n1,n2);
    x = dot(m1,e);
    y = dot(n1,n2);
    
    angle = atan2(x,y);
end

