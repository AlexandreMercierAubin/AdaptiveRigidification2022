function [R] = findVectorRotation(vec)
    %ROTATIONFROMVECTORS takes a 3d vector and generates a rotation matrix
    %that maps to it
    vz = zeros(size(vec));
    vz(3) = 1;
    R = rotationFromVectors(vz,vec);
end

