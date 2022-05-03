classdef Rigidificator3D < handle
    %Rigidificator3D abstract class of the rigidificators which allows the adaptive
    %rigidification of elastic bodies
    properties
       RigidificationThreshold 
       ElastificationThreshold
       RigidificationBendThreshold
       ElastificationBendThreshold
       FrameCount
       Preconditionner = 'ICHOLICT'; %choices; chol or ichol
       Permutation = 'DISSECT';
    end
    
    methods
        function obj = Rigidificator3D()
            %RIGIDIFICATOR Creates a rigidificator, object that decides
            %which part of the mesh are rigid and elastic
            obj.RigidificationThreshold = 1e-3; 
            obj.ElastificationThreshold = 2e-3;
            obj.RigidificationBendThreshold = 1e-3;
            obj.ElastificationBendThreshold = 2e-3;
            obj.FrameCount = 3;
        end
        
        function updateRigidBodies(obj, mesh, rigidTris, settings)
            %UPDATERIGIDBODIES updates the list of rigid bodies given a new
            %list of triangles thatshould be rigid
            ticGraph = tic;
            subgraphedRigid = subgraph(mesh.Graph, rigidTris');
            timeGraph = toc(ticGraph);
            
            
            ticGraph = tic;
            components = conncomp(subgraphedRigid, 'OutputForm', 'cell');
            timeGraph = toc(ticGraph);
            
            
            ticTrans = tic;
            for i = 1 : size(components,2)
                components{i} = rigidTris(components{i})';
            end
            timeTrans = toc(ticTrans);
            
        
            ticUpdateBodies = tic;
            
            % vertices of connected components
            ticVertReshape = tic;
            vertices = cell(1,numel(components));
            for i = 1 : numel(components)
                vertices{i} = reshape(mesh.t(components{i}, :), 1, []);
            end
            timeVertReshape = toc(ticVertReshape);
            
            % find unique values in each components and from each other
            % https://www.mathworks.com/matlabcentral/answers/576337
            ticUnique = tic;
            endIdx = cumsum(cellfun(@numel, vertices));
            startIdx = circshift(endIdx, 1);
            startIdx(1) = 0;
            startIdx = startIdx + 1;
            [C, ixa, ~] = unique([vertices{:}], 'stable');
            timeUnique = toc(ticUnique);
            
            ticRemoveTooSmall = tic;
            vertices = cell(1,size(vertices,2));

            tooSmall = false(1,size(vertices,2));
            
            for i = 1:numel(vertices)
                vertices{i} = C(ixa >= startIdx(i) & ixa <= endIdx(i));
                %No need to check for pinned triangles there
                %vertices{i}(mesh.isTetPinned(vertices{i}) == 1) = [];
                if numel(vertices{i}) < size(mesh.t,2) % 3 in 2D, 4 in 3D
                    tooSmall(i) = true;
                end
            end
            vertices(tooSmall) = [];
            
            timeRemoveTooSmall = toc(ticRemoveTooSmall); 
            
            % get rid of extra bodies we won't need
            ticRemoveExtraBodies = tic;
            mesh.RigidBodies(numel(vertices) + 1:end) = [];
            timeRemoveExtraBodies = toc(ticRemoveExtraBodies);
            
            ticMakeBodies = tic;
            len = numel(mesh.RigidBodies);
            
            if len > 0
                l = min(len, numel(vertices));
                
                % assign vertices of old bodies
                [mesh.RigidBodies(1:l).Indices] = vertices{1:l}; 
                
                % make new required bodies
                mesh.RigidBodies(end + 1:numel(vertices)) = cellfun(@makeBody, vertices(len + 1:numel(vertices)));
            elseif numel(vertices) > 0
                % initialize list of bodies
                mesh.RigidBodies = cellfun(@makeBody, vertices(len + 1:numel(vertices)));
            end
            timeMakeBodies = toc(ticMakeBodies);
            
            % update the bodies
            ticUpdateBodies = tic;
            for body = mesh.RigidBodies
            	 body.updateBody(false);
            end
            timeUpdateBodies = toc(ticUpdateBodies);
            
            if settings.PrintTimings
                disp("Subgraph: " + timeGraph);
                disp("Connected components: " + timeGraph);
                disp("Transpose components: " + timeTrans);
                fprintf('\t\tUpdate bodies: %g\n', toc(ticUpdateBodies));
                disp("vert reshape: " + timeVertReshape );
                disp("unique: " + timeUnique);
                disp("remove too small: " + timeRemoveTooSmall );
                disp("remove extra bodies: " + timeRemoveExtraBodies );
                disp("make bodies: " + timeMakeBodies );
                disp("update bodies: " + timeUpdateBodies );
            end


            ticMesh = tic;
            
            % update the mesh
            mesh.updateRigidState();
            
            if settings.PrintTimings
                fprintf('\t\tUpdate mesh: %g\n', toc(ticMesh));
            end
            
            function body = makeBody(vertices)
                body = RigidBody3D(mesh);
                body.Indices = vertices;
            end
        end
        
        function removeTris(obj, mesh, tris, settings)
            %REMOVETRIS Removes the specified list of triangles from all
            %rigid bodies
            removeTrisTic = tic;
            tic1 = tic;
            
            bodyRefs = mesh.RigidBodies; % list may change
            for x = bodyRefs
                x.removeTriangles(tris, false, settings);
            end
            
            if settings.PrintTimings
                fprintf('\t\tRemove bodies: %g\n', toc(tic1));
            end
            tic2 = tic;
            
            mesh.updateRigidState();
            
            if settings.PrintTimings
                fprintf('\t\tUpdate mesh: %g\n', toc(tic2));
                fprintf('\tRemove tris time: %g\n', toc(removeTrisTic));
            end 
        end
        
    end
    
    methods (Abstract = true)      
        checkForElastification( obj, mesh, cache, frame, h, settings );
    end
end

