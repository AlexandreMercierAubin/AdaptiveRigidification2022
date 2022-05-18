function singleUseFigures = plotMesh( h, mesh2d, settings, singleUseFigures, cache, integrator, rigidificator)
    % plotMesh plots given edges, mesh, and pinned vertices (assumes hold on)
   
    if nargin < 3
        settings = SimulationSettings(); 
    end
    
    pr = reshape(mesh2d.p, 2, mesh2d.N)';
        
    edgeCol = iif(settings.DrawEdges, 2*[2 1 34] / 255, 'none');
    faceMaterials = mesh2d.materials(mesh2d.materialIndex);
    faceColor = [faceMaterials.color];
    faceColor = reshape(faceColor,3,[])';
                
    % Add a red tint to rigid parts
%     tintCol = [0.8 0.5 0.2];  % orange
    tintCol = [0.8 0.0 0.0];  % red
    if settings.plotTriImplicitRigidificationElastification
        tintColImplicit = [0.4,0.8,0.4]; %when -1
        tintColRemoved = [0.8,0.4,0.8]; %when 1
    end
    tintFactor = 0.5;
    if isa(mesh2d, 'AdaptiveMesh')
        for j = 1:numel(mesh2d.RigidBodies)
            inds = mesh2d.RigidBodies(j).TriInds';
            faceColor(inds, :) = faceColor(inds, :) + (tintCol - faceColor(inds, :)) * tintFactor ;
        end
        
        if settings.plotTriImplicitRigidificationElastification
            indsPos = mesh2d.RigidificationDifference == 1;
            indsNeg = mesh2d.RigidificationDifference == -1;
            faceColor(indsPos, :) = faceColor(indsPos, :) + (tintColRemoved - faceColor(indsPos, :)) * tintFactor ;
            faceColor(indsNeg, :) = faceColor(indsNeg, :) + (tintColImplicit - faceColor(indsNeg, :)) * tintFactor ;
        end
    end
    
    
    if mesh2d.RenderPatch == 0
        mesh2d.RenderPatch = patch('vertices', pr, 'faces', mesh2d.t, 'edgecol', edgeCol, 'facecol', 'flat', 'FaceVertexCData', faceColor, 'FaceAlpha', .5, 'EdgeAlpha', .03 );    
    else
        mesh2d.RenderPatch.Vertices = pr;
        mesh2d.RenderPatch.FaceVertexCData = faceColor;
    end
        
    %plot(pr(mesh.pinnedInds, 1), pr(mesh.pinnedInds, 2), 'r.');

    % 
    if ( numel(mesh2d.BoundaryLines) == 0 )
        mesh2d.BoundaryLines = cell(1,numel(mesh2d.boundaryEdges));
        for j=1:numel(mesh2d.boundaryEdges)
            edges = mesh2d.boundaryEdges{j};
            boundaryx = [ pr(edges(:,1),1); pr(edges(1,1),1) ];
            boundaryy = [ pr(edges(:,1),2); pr(edges(1,1),2) ];
            mesh2d.BoundaryLines{j} = plot( boundaryx, boundaryy, 'k-' );
            axis off;
        end
    else
        for j=1:numel(mesh2d.boundaryEdges)
            edges = mesh2d.boundaryEdges{j};
            boundaryx = [ pr(edges(:,1),1); pr(edges(1,1),1) ];
            boundaryy = [ pr(edges(:,1),2); pr(edges(1,1),2) ];    
            mesh2d.BoundaryLines{j}.XData = boundaryx;
            mesh2d.BoundaryLines{j}.YData = boundaryy;
        end
    end
        
    if settings.DrawElasticForces && numel(cache.elasticForces) > 0
        color = [0.7,0.1,0.1];
        inds = find(~mesh2d.pinned);
        for j=1:numel(inds)
            i = inds(j)*2;
            posx = [mesh2d.p(i-1);mesh2d.p(i-1) + cache.elasticForces(i-1)];
            posy = [mesh2d.p(i);mesh2d.p(i) + cache.elasticForces(i)];
            t = line(posx,posy, 'color',color');
            singleUseFigures = [singleUseFigures,t'];
        end
    end
    
    constTimes = 10;
    if settings.DrawDv
        color = [1,0,0];
        for i = 1:2:numel(cache.oldDv)
            posx = [mesh2d.p(i);mesh2d.p(i) + constTimes*cache.oldDv(i)];
            posy = [mesh2d.p(i+1);mesh2d.p(i+1) + constTimes*cache.oldDv(i+1)];
            t = line(posx,posy, 'color',color','LineWidth',3);
            singleUseFigures = [singleUseFigures,t'];
        end
    end
    
    if settings.DrawApproxDv
        color = [0,0.7,0]; % slightly darker green
        for i = 1:2:numel(cache.ApproximatedDeltaV)
            posx = [mesh2d.p(i);mesh2d.p(i) + constTimes*cache.ApproximatedDeltaV(i)];
            posy = [mesh2d.p(i+1);mesh2d.p(i+1) + constTimes*cache.ApproximatedDeltaV(i+1)];
            t = line(posx,posy, 'color',color','LineWidth',1);
            singleUseFigures = [singleUseFigures,t'];
        end
    end
    % Everything below is deprecated, and should be made in a way that do
    % not add elements to the plot every time. In practice, it is still
    % usable, although you will need to clear the plot everytime. I don't 
    % recommend doing that because it is very slow, but its OK for
    % debugging.
    
    if settings.DrawContactFrames
        n = numel(cache.cInfo);
        for i=1:n
            contactInfo = cache.cInfo(i);

            p = contactInfo.point;
            pos = [p; p + contactInfo.normal];
            
            color = [1,0,0];
            n = line(pos(:,1),pos(:,2), 'color',color');
            pos = [p; p + contactInfo.tangent];
            t = line(pos(:,1),pos(:,2), 'color',color');
            
            singleUseFigures = [singleUseFigures,n',t'];
            if ~isnan(contactInfo.pointAlpha)
                color = [0,1,0];
                p = contactInfo.pointAlpha;
                pos = [p; p + contactInfo.normalAlpha];
                n = line(pos(:,1),pos(:,2), 'color',color');
                pos = [p; p + contactInfo.tangentAlpha];
                t = line(pos(:,1),pos(:,2), 'color',color');
                singleUseFigures = [singleUseFigures,n',t'];
                
                color = [1,0,1];
                p = contactInfo.pointAlphaInv;
                pos = [p; p + contactInfo.normalAlphaInv];
                n = line(pos(:,1),pos(:,2), 'color',color);
                pos = [p; p + contactInfo.tangentAlphaInv];
                t = line(pos(:,1),pos(:,2), 'color',color');
                singleUseFigures = [singleUseFigures,n',t'];
            end
            
            
        end
    end
    if settings.DrawForces
        hh = h * 10;

        for i = 1:mesh2d.N
            % what if we try drawing all the forces?
            if ismember(i, mesh2d.pinnedInds) % || isa(mesh, 'AdaptiveMesh') && ( ismember(i, mesh.ElasticInds) == 0 )
                continue;
            end
            mesh2d.f = zeros(mesh2d.N*2,1);
            mesh2d.f(2:2:end) = mesh2d.mass(2:2:end) * integrator.Gravity;
            if numel(cache.elasticForces)>0
                mesh2d.f = mesh2d.f + cache.elasticForces + integrator.CustomForce;
            end
            
            f = fastQuiver( [mesh2d.p(i * 2 - 1, 1), mesh2d.p(i * 2, 1)], [hh * mesh2d.f(i * 2 - 1, 1), hh * mesh2d.f(i * 2, 1)],  [252 158 79] / 255 /2);
            singleUseFigures = [singleUseFigures,f'];
        end
        
        if isa(mesh2d, 'AdaptiveMesh')
            hh = h * 1;
            for i = 1:numel(mesh2d.RigidBodies)
                body = mesh2d.RigidBodies(i);
                f = fastQuiver( [ body.Position(1), body.Position(2)], [hh * body.Force(1), hh * body.Force(2)], [252 0 0] / 255/2 );                
                %cf = circular_arrow(gcf, 0.5, body.Position, rad2deg(body.Angle), min(330, h * abs(rad2deg(body.Torque))), -sign(body.Torque), [252 0 0] / 255, 5);
                singleUseFigures = [singleUseFigures,f'];
            end
        end
        
    end
    
    if settings.DrawVelocities
        for i = 1:mesh2d.N
            if ismember(i, mesh2d.pinnedInds) || isa(mesh2d, 'AdaptiveMesh') && ismember(i, mesh2d.ElasticInds) == 0
                f = fastQuiver( [ mesh2d.p(i * 2 - 1, 1), mesh2d.p(i * 2, 1)], [mesh2d.v(i * 2 - 1, 1), mesh2d.v(i * 2, 1)], [252 0 0] / 255/2 );                
                singleUseFigures = [singleUseFigures,f'];
                continue;
            else
                f = fastQuiver( [ mesh2d.p(i * 2 - 1, 1), mesh2d.p(i * 2, 1)], [mesh2d.v(i * 2 - 1, 1), mesh2d.v(i * 2, 1)], [84 140 47] / 255/2 );                
                singleUseFigures = [singleUseFigures,f'];
            end
        end
        
         if isa(mesh2d, 'AdaptiveMesh')
             for i = 1:numel(mesh2d.RigidBodies)
                 body = mesh2d.RigidBodies(i);
                 f = fastQuiver([body.Position(1), body.Position(2)], [body.Velocity(1), body.Velocity(2)], [84 140 47] / 255);
                 cf = circular_arrow(gcf, 0.5, body.Position, rad2deg(body.Angle), min(330, abs(rad2deg(body.AngularVelocity))), -sign(body.AngularVelocity), [84 140 47] / 255, 5);
                 singleUseFigures = [singleUseFigures,f',cf];
             end
         end

    end
    
    if settings.DrawLagrangeMultipliers
        for i = 1:size(mesh2d.t, 1)
            
            lambdaX1 = mesh2d.lagrangeMults(i * 4 - 3);
            lambdaY1 = mesh2d.lagrangeMults(i * 4 - 2);
            len1 = sqrt(lambdaX1 * lambdaX1 + lambdaY1 * lambdaY1);

            lambdaX2 = mesh2d.lagrangeMults(i * 4 - 1);
            lambdaY2 = mesh2d.lagrangeMults(i * 4);
            len2 = sqrt(lambdaX2 * lambdaX2 + lambdaY2 * lambdaY2);
            
            pX = 0;
            pY = 0;
            
            for j = 1:3
                index = mesh2d.t(i, j);
                
                pX = pX + mesh2d.p(index * 2 - 1, 1) / 3;
                pY = pY + mesh2d.p(index * 2, 1) / 3;
            end

            v1 = cheap_pseudo_sigmoid(abs(len1));
            v2 = cheap_pseudo_sigmoid(abs(len2));
            
            f = fastQuiver([pX, pY], [lambdaX1, lambdaY1], [v1 0 0]);
            g = fastQuiver([pX, pY], [lambdaX2, lambdaY2], [v2 0 0]);
            singleUseFigures = [singleUseFigures,f',g'];
        end
    end
    
    if isa(mesh2d, 'AdaptiveMesh')
        if settings.DrawRigidFrames
            for i = 1:numel(mesh2d.RigidBodies)
                body = mesh2d.RigidBodies(i);
                f = fastQuiver(body.Position', [cos(body.Angle), sin(body.Angle)], [1 0 1]);
                g = fastQuiver(body.Position', [cos(body.Angle + pi / 2), sin(body.Angle + pi / 2)], [1 0 1]); 
                singleUseFigures = [singleUseFigures,f',g'];
            end
        end
        
        if settings.DrawEDots && ~isempty(mesh2d.RigidificationValues)           
            Fa = mesh2d.B * mesh2d.p;
            Fb = mesh2d.B * cache.oldp; % not so efficient, these would have been previously computed
            
            EDotNorms = mexEdiffNorm2D( Fa, Fb, h );
            
            for t = 1:size(mesh2d.t)
                i = mesh2d.t(t, 1);
                j = mesh2d.t(t, 2);
                k = mesh2d.t(t, 3);
                
                center(:,t) = mesh2d.p(i * 2 - 1:i * 2) + mesh2d.p(j * 2 - 1:j * 2) + mesh2d.p(k * 2 - 1:k * 2);
                center(:,t) = center(:,t) / 3;
            end
            color = -log(EDotNorms); %norm(mesh2d.B(:,8*2-1:8*2)*mesh2d.p0(8*2-1:8*2),'fro')
            color = color./8;
            color = max(min(color,1),0);
            f = scatter(center(1,:), center(2,:), 5, [color,color,color], 'filled');
            singleUseFigures = [singleUseFigures,f];
            
        end
        
        if ~settings.DrawApproxEDots
            if isfield( mesh2d.reusableFigs, 'approxEDots' ) && isvalid( mesh2d.reusableFigs.approxEDots )
                delete( mesh2d.reusableFigs.approxEDots );
            end
        elseif settings.DrawApproxEDots && ~isempty(cache.edotnorms)
            vel = mesh2d.v + cache.ApproximatedDeltaV; 
            F = mesh2d.B * (mesh2d.p + h*vel+mesh2d.p)/2;
            Fdot = mesh2d.B * vel;

            EDotNorms = mexEdotNorm2D( F, Fdot );

            p = reshape( mesh2d.p, 2, [] );
            center = (p(:,mesh2d.t(:,1)) + p(:,mesh2d.t(:,2)) + p(:,mesh2d.t(:,3))) / 3;

            color = max(min(1-log(EDotNorms),1),0); %norm(mesh2d.B(:,8*2-1:8*2)*mesh2d.p0(8*2-1:8*2),'fro')
            elastify = EDotNorms <= rigidificator.ElastificationThreshold; %modulate color of tagged for elastification
            
            %RUN this code in the console command line with a breakpoint here 
            %if you want to plot the error with respect to identity
            
%             error = zeros(size(mesh2d.t,1),1);
%             sumF = zeros(2,2);
%             for t = 1:size(mesh2d.t)
%                 Brows = [t*4-3,t*4-2,t*4-1,t*4];
%                 verts = mesh2d.t(t,:);
%                 Bcols = [verts(1)*2-1,verts(1)*2,verts(2)*2-1,verts(2)*2,verts(3)*2-1,verts(3)*2];
%                 localB = mesh2d.B(Brows,Bcols);
%                 localp0 = mesh2d.p0(Bcols);
%                 localF = reshape(localB*localp0,2,2);
%                 sumF = sumF + localF;
%                 respectI = localF-eye(2);
%                 error(t) = norm(respectI, 'fro');
%             end
% %             figure;plot(error);
%             meanF = sumF./size(mesh2d.t,1);
%             for t = 1:size(mesh2d.t)
%                 Brows = [t*4-3,t*4-2,t*4-1,t*4];
%                 verts = mesh2d.t(t,:);
%                 Bcols = [verts(1)*2-1,verts(1)*2,verts(2)*2-1,verts(2)*2,verts(3)*2-1,verts(3)*2];
%                 localB = mesh2d.B(Brows,Bcols);
%                 localp0 = mesh2d.p0(Bcols);
%                 localF = reshape(localB*localp0,2,2);
%                 respectMean = localF-meanF;
%                 errorMean(t) = norm(respectMean, 'fro');
%             end
%             figure;plot(errorMean);

            if isfield( mesh2d.reusableFigs, 'approxEDots' ) && isvalid( mesh2d.reusableFigs.approxEDots )
                mesh2d.reusableFigs.approxEDots.XData = center(1,:);
                mesh2d.reusableFigs.approxEDots.YData = center(2,:);
                mesh2d.reusableFigs.approxEDots.CData = [color,color.*elastify,color.*elastify];
            else
                mesh2d.reusableFigs.approxEDots = scatter(center(1,:), center(2,:), 10, [color,color.*elastify,color.*elastify], 'filled');
            end

        end
        
        if settings.DrawEigensOfE
            for t = 1:size(mesh2d.t)
                i = mesh2d.t(t, 1);
                j = mesh2d.t(t, 2);
                k = mesh2d.t(t, 3);
                
                p_i = mesh2d.p(i * 2 - 1:i * 2);
                p_j = mesh2d.p(j * 2 - 1:j * 2);
                p_k = mesh2d.p(k * 2 - 1:k * 2);
                
                Ds = [p_i - p_k, p_j - p_k];
                F = Ds * mesh2d.el(t).Bm;
                
                E = 0.5 * (F' * F - eye(2));
                
                center = mesh2d.p(i * 2 - 1:i * 2) + mesh2d.p(j * 2 - 1:j * 2) + mesh2d.p(k * 2 - 1:k * 2);
                center = center / 3;
                
                [V, D] = eig(E);
                eigenVals = diag(D);
                
                for i = 1:numel(eigenVals)
                    val = eigenVals(i);
                    vec = V(:, i);
                    f = fastQuiver(center, vec * val, [1 0 0]);
                    singleUseFigures = [singleUseFigures,f'];
                end
            end
        end
    end
end