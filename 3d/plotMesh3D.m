function [] = plotMesh3D( mesh3D, cache, settings)
% plotMesh plots given edges, mesh, and pinned vertices (assumes hold on)

pr = reshape(mesh3D.p, 3, mesh3D.N)';
plotFacesOnly = true;
edgeCol = [2 1 34] / 255;
if plotFacesOnly && isempty(mesh3D.SurfaceWeights)
    materialList = mesh3D.materials(mesh3D.faceMaterialIndex);
else
    materialList = mesh3D.materials(mesh3D.materialIndex);
end

elementColor = zeros(numel(materialList),3);
for i = 1:numel(materialList)
    elementColor(i,:) = materialList(i).color;
end
faceColor = elementColor;

faceAlpha = 1;

% Add a red tint to rigid parts
tintCol = [0.8 0.5 0.2];
tintFactor = 0.75;

activatedInds = false(size(elementColor,1),1);
if isa(mesh3D, 'AdaptiveMesh3D')
    if plotFacesOnly && isempty(mesh3D.SurfaceWeights)
        for j = 1:numel(mesh3D.RigidBodies)
            inds = mesh3D.RigidBodies(j).TetInds';
            faceInds = find(ismember(mesh3D.facesTriList,inds));
            activatedInds(faceInds) = true;
        end
    else
        activatedInds = ~mesh3D.isTetElastic;
    end
    faceColor(activatedInds, :) = elementColor(activatedInds, :) + (tintCol - elementColor(activatedInds, :)) * tintFactor ;
    
    
    if settings.DrawRigidDOFs
        faceAlpha = 0.5;
        rigidDofs = [mesh3D.RigidBodies.Indices];
        
        for i = 1 : mesh3D.N
            c = [mesh3D.p(3*i-2),mesh3D.p(3*i-1),mesh3D.p(3*i)];
            r = 0.05;

            [x,y,z] = sphere;
            x = x * r;
            y = y * r;
            z = z * r;

            xdata = c(1)+x;
            ydata = c(2)+y;
            zdata = c(3)+z;
            
            cdata = ones(size(zdata));
            cdata(1) = 0;
            
            hold on;
            if numel(mesh3D.DOFsSpheres) < mesh3D.N
                plotHandle = surf(xdata,ydata,zdata,cdata,'EdgeColor','none','FaceColor',[0 0 0]);
                mesh3D.DOFsSpheres = [mesh3D.DOFsSpheres,plotHandle];
            else
                mesh3D.DOFsSpheres(i).XData = xdata;
                mesh3D.DOFsSpheres(i).YData = ydata;
                mesh3D.DOFsSpheres(i).ZData = zdata;
                set(mesh3D.DOFsSpheres(i),'FaceColor',[0 0 0]);
            end
        end
        
        for i = 1:numel(rigidDofs)
            set(mesh3D.DOFsSpheres(rigidDofs(i)),'FaceColor',[1 0 0]);
        end
    end
end
if ~isempty(mesh3D.SurfaceWeights)
    V = reshape(mesh3D.p',3,[])';
    nV = mesh3D.SurfaceWeights * V;
    surfaceFaceColor = faceColor(mesh3D.Surface2Tet,:);
    if mesh3D.RenderPatch == 0
        mesh3D.RenderPatch = patch('vertices', nV, 'faces', mesh3D.SurfaceFaces, 'edgecol', [0.2, 0.2, 0.2],'facecol', 'flat', 'FaceVertexCData', surfaceFaceColor , 'FaceAlpha', faceAlpha, 'EdgeColor','none'); 
    else
        mesh3D.RenderPatch.Vertices = nV;
        mesh3D.RenderPatch.FaceVertexCData = surfaceFaceColor;
    end
else
    if mesh3D.RenderPatch == 0
        if plotFacesOnly
            mesh3D.RenderPatch = patch('vertices', pr, 'faces', mesh3D.faces, 'edgecol', [0.2, 0.2, 0.2], 'facecol', 'flat', 'FaceVertexCData', faceColor , 'FaceAlpha', faceAlpha, 'EdgeColor','none'); 
        else
            mesh3D.RenderPatch = patch('vertices', pr, 'faces', mesh3D.t, 'edgecol', [0.2, 0.2, 0.2], 'facecol', 'flat', 'FaceVertexCData', faceColor , 'FaceAlpha', .05, 'EdgeColor', 'none'); 
        end
    else
        mesh3D.RenderPatch.Vertices = pr;
        mesh3D.RenderPatch.FaceVertexCData = faceColor;
    end
end

if settings.DrawVelocities
    for i = 1:mesh3D.N
        xdata = [mesh3D.p(i*3-2);mesh3D.p(i*3-2)+mesh3D.v(i*3-2)];
        ydata = [mesh3D.p(i*3-1);mesh3D.p(i*3-1)+mesh3D.v(i*3-1)];
        zdata = [mesh3D.p(i*3);mesh3D.p(i*3)+mesh3D.v(i*3)];
        if numel(mesh3D.lineHandles) < i
            mesh3D.lineHandles = [mesh3D.lineHandles;line(xdata,ydata,zdata, 'Color',[100/255,100/255,255/255])];
        else
            mesh3D.lineHandles(i).XData = xdata;
            mesh3D.lineHandles(i).YData = ydata;
            mesh3D.lineHandles(i).ZData = zdata;
        end
    end
elseif settings.DrawForces
    for i = 1:mesh3D.N
        xdata = [mesh3D.p(i*3-2);mesh3D.p(i*3-2)+mesh3D.f(i*3-2)];
        ydata = [mesh3D.p(i*3-1);mesh3D.p(i*3-1)+mesh3D.f(i*3-1)];
        zdata = [mesh3D.p(i*3);mesh3D.p(i*3)+mesh3D.f(i*3)];
        if numel(mesh3D.lineHandles) < i
            mesh3D.lineHandles = [mesh3D.lineHandles;line(xdata,ydata,zdata, 'Color',[100/255,100/255,255/255])];
        else
            mesh3D.lineHandles(i).XData = xdata;
            mesh3D.lineHandles(i).YData = ydata;
            mesh3D.lineHandles(i).ZData = zdata;
        end
    end
elseif (settings.DrawLambdas && ~isempty(cache.cInfo))    
    cInfos = cache.cInfo;
    for j = 1:numel(cInfos)
        cp = cInfos(j).point;
        if settings.DrawLambdas
            shorten = settings.DrawLambdasScale;
            cn = cInfos(j).normal;
            ct = cInfos(j).tangent;
            ct2 = cInfos(j).tangent2;
            cl = cache.prevLambdas(j*3-2:j*3);
            cf = (cn*cl(1) + ct*cl(2) + ct2*cl(3)) * shorten;
            xdata = [cp(1); cp(1)-cf(1)];
            ydata = [cp(2); cp(2)-cf(2)];
            zdata = [cp(3); cp(3)-cf(3)];
            if numel(mesh3D.lineHandles) < j
                mesh3D.lineHandles = [mesh3D.lineHandles;line(xdata,ydata,zdata, 'Color',[100/255,100/255,255/255])];
            else
                mesh3D.lineHandles(j).XData = xdata;
                mesh3D.lineHandles(j).YData = ydata;
                mesh3D.lineHandles(j).ZData = zdata;
            end
        end
    end
    %hides the extra lambdas from contacts that are no longer there
    for j = numel(cInfos):numel(mesh3D.lineHandles)
            xdata = [0; 0];
            ydata = [0; 0];
            zdata = [0; 0];
            mesh3D.lineHandles(j).XData = xdata;
            mesh3D.lineHandles(j).YData = ydata;
            mesh3D.lineHandles(j).ZData = zdata;
    end
end
    
end


