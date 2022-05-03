function applyAcceleration(mesh, acc)
    % Probably doesn't matter if we apply forces to pinned particles...
    % they'll not be updated or moved!  But be careful!
    mesh.f = repmat( acc, mesh.N, 1) .* mesh.mass .* (1-reshape( repmat( mesh.pinned', 3, 1 ), mesh.N*3, 1 ) );
end

