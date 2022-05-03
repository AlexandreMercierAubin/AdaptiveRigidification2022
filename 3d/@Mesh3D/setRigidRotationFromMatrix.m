function  setRigidRotationFromMatrix( mesh3d, R, forceCenterOfMassRotation )
%RIGIDTRANSFORM Transforms mesh points by given rotation and translation.
%   mesh.rigidTransform( degrees, translation, scale )
%   degrees should be a 3 component vector axis x, y ,z
    if nargin < 4
        forceCenterOfMassRotation = false;
    end
        
    if ~forceCenterOfMassRotation
        p = [ mesh3d.p(1:3:end), mesh3d.p(2:3:end), mesh3d.p(3:3:end)]';
        mesh3d.p(:) = R * p;
    else
        masses = mesh3d.mass(1:3:end);
        mass = sum(masses);
        centerOfMass = [dot(mesh3d.p(1:3:end), masses); dot(mesh3d.p(2:3:end), masses); dot(mesh3d.p(3:3:end), masses)] / mass;
        p = reshape( mesh3d.p', 3, [] )' - centerOfMass';
        newp = R * p';
        newp = newp' + centerOfMass';
        mesh3d.p = reshape(newp',1,[])';
    end

end

