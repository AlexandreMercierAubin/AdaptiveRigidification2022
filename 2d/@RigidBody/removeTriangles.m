function removeTriangles( obj, tris, updateMesh, settings )
    %REMOVETRIANGLES Removes triangles in the body
    %   Given a list of triangles, removes the triangles and
    %   disconnect the body into multiple rigid bodies if necessary. This
    %   function may delete the body and alter the number of body in the mesh
    %   dramatically. Ignores triangles that are not part of the body.
    
    if nargin < 3
        updateMesh = 1;
    end

    % find new vertices
    %a logical setdiff
    triangleIdRange = [1:obj.Mesh.N]';
    setDiffIds = false(numel(triangleIdRange),1);
    setDiffIds(obj.Indices) = true;
    setDiffIds(obj.Mesh.t(tris, :)) = false;
    vertices = triangleIdRange(setDiffIds);
%     vertices = setdiff(obj.Indices, [obj.Mesh.t(tris, :)]);

    if numel(vertices) == numel(obj.Indices)
        return;
    end

    % find new tris
    if numel(vertices) < 3
        tris = [];
    else
        % the triangle index must be found 3 times, once for each adjacent
        % vertex, if it is to remain part of a rigid body
        tris = find(accumarray([obj.Mesh.TrianglesPerParticle{vertices}]', 1) == 3)';
    end
    
    if isempty(tris)
        ticDelete = tic;
        
        obj.Mesh.RigidBodies(obj.Mesh.RigidBodies == obj) = [];
        
        if settings.PrintTimings
            fprintf('\t\t\tDelete empty body: %g\n', toc(ticDelete));
        end
        return;
    end

    ticGraph = tic;
    
    % matlab's builtin is thankfully very fast
    components = conncomp(subgraph(obj.Mesh.Graph, tris), 'OutputForm', 'cell');
    components = cellfun(@(x) tris(x), components, 'UniformOutput', false);
    
    if settings.PrintTimings
        fprintf('\t\t\tFind connected components: %g\n', toc(ticGraph));
    end
    ticRemove2 = tic;
    
    % get the vertices of each components
    vertices = cellfun(@(x) reshape(obj.Mesh.t(x, :), 1, []), components, 'UniformOutput', false);
    
    % this whole thing just removes dupplicates from each of the vertices
    % cells and from each other.
    % more info: https://www.mathworks.com/matlabcentral/answers/576337
    endIdx = cumsum(cellfun(@numel, vertices));
    startIdx = circshift(endIdx, 1);
    startIdx(1) = 0;
    startIdx = startIdx + 1;
    [C, ixa, ~] = unique([vertices{:}], 'stable');

    for i = 1:numel(vertices)
        vertices{i} = C(ixa >= startIdx(i) & ixa <= endIdx(i));
    end
    tooSmall = cellfun(@(x) numel(x) < 3, vertices);  %can this ever happen?
    vertices(tooSmall) = [];
    
    % update indices
    obj.Indices = vertices{1};
    
    if numel(vertices) > 1
        % make more bodies for the rest of connected components
        newBodies = cellfun(@makeBody, vertices(2:end));
        obj.Mesh.RigidBodies = [obj.Mesh.RigidBodies, newBodies];
        
        for i = newBodies
            i.updateBody(false); % update them, not updating the mesh
        end
    end
    
    % update this body, only updating the mesh if necessary
    obj.updateBody(updateMesh);
    if settings.PrintTimings
        fprintf('\t\t\tRemove triangles: %g\n', toc(ticRemove2));
    end
    
    function body = makeBody(vertices)
        body = RigidBody(obj.Mesh);
        body.Indices = vertices;
    end
end