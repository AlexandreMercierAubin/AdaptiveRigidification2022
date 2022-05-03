function outputSceneRenderFile(mesh3D, contactFinders, settings, k, frame)
   V = reshape(mesh3D.p, 3, [])';
   F = mesh3D.faces;
   
   %surface mesh
   if ~isempty(mesh3D.SurfaceWeights)
        nV = mesh3D.SurfaceWeights * V;
        if isa(mesh3D, "AdaptiveMesh3D")
            rigidSurfaceInds = ~mesh3D.isTetElastic(mesh3D.Surface2Tet);
        else
            rigidSurfaceInds = [];
        end
        VColor = mesh3D.surfaceVertexBaseColor;
        tintCol = [0.8 0.5 0.2];
        tintFactor = 0.75;
        VColor(rigidSurfaceInds,:) = VColor(rigidSurfaceInds,:) + tintCol - VColor(rigidSurfaceInds,:) * tintFactor;
        VColor = floor(VColor * 255);
        filenamePly = [settings.OBJDir 'frame_surface_' num2str(k,'%03.f') '_' num2str(frame/(settings.PlotSkip+1),'%03.f'), '.ply'];
        writePLY(filenamePly, nV, mesh3D.SurfaceFaces,'ascii', VColor);
%                         writeOBJ([settings.OBJDir 'frame_surface_' num2str(k,'%03.f') '_' num2str(frame/(settings.PlotSkip+1),'%03.f'), '.obj'], nV, mesh3D.SurfaceFaces);
   end
   
   if settings.ExportNormalsOBJs
       N = normals(V,F, 'Stable', false, 'UseSVD', true);
       for i = 1 : size(N,1)
           normN = norm(N(i,:));
           N(i,:) = N(i,:)/normN;
       end
       writeOBJ([settings.OBJDir 'frame_' num2str(k,'%03.f') '_' num2str(frame/(settings.PlotSkip+1),'%03.f'), '.obj'], V, F, [], [],N);
   else
       writeOBJ([settings.OBJDir 'frame_' num2str(k,'%03.f') '_' num2str(frame/(settings.PlotSkip+1),'%03.f'), '.obj'], V, F);
   end

   %write rigid body chunk
   if isa(mesh3D, 'AdaptiveMesh3D')
       inds = [];
       for m = 1:numel(mesh3D.RigidBodies)
           inds = [inds; mesh3D.RigidBodies(m).TetInds'];
       end
       inds = unique(inds);
       %faceInds = find(ismember(mesh3D(j).facesTriList,inds));
       faceInds = boundary_faces(mesh3D.t(inds,:));
       if isempty(faceInds) == 1
           faceInds = [1 1 1];
       end

       %rigid part
       V = reshape(mesh3D.p, 3, [])';
       F = faceInds;
       if settings.ExportNormalsOBJs
           N = normals(V,F, 'Stable', false, 'UseSVD', true);
           for i = 1 : size(N,1)
               normN = norm(N(i,:));
               N(i,:) = N(i,:)/normN;
           end
           writeOBJ([settings.OBJDir 'frame_rigid_' num2str(k,'%03.f') '_' num2str(frame/(settings.PlotSkip+1),'%03.f'), '.obj'], V, F, [],[], N);
       else
           writeOBJ([settings.OBJDir 'frame_rigid_' num2str(k,'%03.f') '_' num2str(frame/(settings.PlotSkip+1),'%03.f'), '.obj'], V, F);
       end
   end

   %write the contact detector renders
   contactV = [];
   contactF = [];
   for j = 1:numel(contactFinders)
        [V,F] = contactFinders{j}.getObjPositionFaces(frame);
        N = size(contactV,1);
        contactV = [contactV;V];
        contactF = [contactF;N+F];
   end
   if numel(contactFinders) > 0
       writeOBJ([settings.OBJDir 'frame_contactRender_' num2str(k,'%03.f') '_' num2str(frame/(settings.PlotSkip+1),'%03.f'), '.obj'], contactV, contactF);
   end 
end

