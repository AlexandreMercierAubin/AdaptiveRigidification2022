function R = expRodrigues( omega, h )
    %EXPRODRIGUES returns a rotation matrix exp(h * omegahat);
    % omega     3D angular velocity
    % h         time step
    
    mag = norm( omega );
    w = omega / mag;
    
    if mag < 1e-6
        R = eye(3);
        return;
    end
    
    theta = h * mag;
    c = cos(theta);
    s = sin(theta);

    c1 = 1 - c;
    R = zeros(3,3);

    % Hard to read this and know it is correct... 
    R(1,1) =  c        + w(1) * w(1) * c1;
    R(2,1) =  w(3) * s + w(1) * w(2) * c1;
    R(3,1) = -w(2) * s + w(1) * w(3) * c1;

    R(1,2) = -w(3) * s + w(1) * w(2) * c1;
    R(2,2) =  c        + w(2) * w(2) * c1;
    R(3,2) =  w(1) * s + w(2) * w(3) * c1;

    R(1,3) =  w(2) * s + w(1) * w(3) * c1;
    R(2,3) = -w(1) * s + w(2) * w(3) * c1;
    R(3,3) =  c        + w(3) * w(3) * c1;
end