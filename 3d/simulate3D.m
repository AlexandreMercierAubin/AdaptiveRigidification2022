% PHYSICS PARAMS
function td = simulate3D(meshes, h, contactFinders, integrators, rigidificators, settings, energyModel, animationScripters)
    
    % ensures a valid number of input arguments
    if nargin < 8
        animationScripters = NullAnimationScript();
    end

    %turn the inputs into cell arrays
    if ~isa(animationScripters, 'cell')
        animationScripters = { animationScripters };
    end
    
    if ~isa(meshes, 'cell')
        meshes = { meshes };
    end
    
    if ~isa(integrators, 'cell')
        integrators = { integrators };
    end
    
    if ~isa(rigidificators, 'cell')
        rigidificators = { rigidificators };
    end
    
    if ~isa(contactFinders, 'cell')
        contactFinders = { contactFinders };
    end

    [comparisons, initialMeshes, meshes, integrators, rigidificators] = setupComparisons(meshes, integrators, rigidificators);

    %UI variables TODO: make a UI class so this is clean
    mouseDown = 0;
    mousePoint = [nan, nan, nan]; % Current mouse point
    pullDirection = [nan nan nan];
    mouseGrabMeshId = nan;  % Grabbed mesh
    mouseLine = 0;

    %create timing data
    td = cell( comparisons, 1 );
    for k = 1:comparisons
        td{k} = TimingData();
    end
    
    [mainFig,axesList, initialCamera] = setupWindow(settings,meshes);
    set(mainFig, 'WindowKeyPressFcn', @onKeyPressed);
    set(mainFig, 'WindowButtonUpFcn', @onMouseUp);
    set(mainFig, 'WindowButtonDownFcn', @onMouseDown);
    set(mainFig, 'WindowButtonMotionFcn',@onMouseDrag);

    caches = setupCache(settings, meshes, rigidificators, comparisons, energyModel, h);
    

    if settings.MakeVideo == 1
        v = VideoWriter(strcat('out/', convertStringsToChars(settings.SceneName),'.mp4'),'MPEG-4');
        open(v);
    end

    running = settings.RunRightAway;
    stepOnce = 0;
    notDone = true;

    frame = 0;
    elapsed = 0;
    skip = settings.PlotSkip + 1;
    % tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'normal'); 
    tStart = cputime;
    
    while ishghandle(mainFig) && notDone

        if frame >= settings.FramesToRecord && settings.FramesToRecord ~= -1
            notDone = false;
        end
        %(settings.FramesToRecord == -1 || settings.FramesToRecord > frame)
        % If we are not running, or don't have a step request, then we
        % really shouldn't bother doing all the redrawing!  Only collect
        % video frames for when the simulation is running!
        if ~running
            if stepOnce
                disp('Stepping once');
                stepOnce = 0; 
            else
                pause(0.2);
                continue;
            end
        end 


        set(0, 'CurrentFigure', mainFig)

        

        %% -------- BEGIN DRAWING -----    
        frame = frame + 1;
        if settings.printFrameNumber
            frame
        end
        if mod(frame, skip) == 0
            hold on;
            
            for k = 1:comparisons
                mesh3D = meshes{k};
                rigidificator = rigidificators{k};
                integrator = integrators{k};
                cache = caches{k};
                
                for j = 1:numel(contactFinders)
                    contactFinders{j}.render(frame);
                end
                
                plotMesh3D( mesh3D, cache, settings );  
                %title(sprintf('t = %f', elapsed ));
                % the title is not visible with the given camera... the elapsed
                % time could perhaps be drawn onto the figure in another way

                if settings.WriteOBJs == 1
                    outputSceneRenderFile(mesh3D, contactFinders, settings, k, frame);
                end
            end
            if settings.MakeVideo == 1
                F = getframe(gcf);
                writeVideo(v,F);
            else
                drawnow;
            end
        end
        for k = 1:comparisons
            mesh3D = meshes{k};
            rigidificator = rigidificators{k};
            integrator = integrators{k};
            cache = caches{k};

            %when specified, focus the camera on a vertex
            if ~isempty(settings.FocusOnMeshNode)
                xPos = mesh3D.p(settings.FocusOnMeshNode(2)*3-2);
                yPos = mesh3D.p(settings.FocusOnMeshNode(2)*3-1);
                zPos = mesh3D.p(settings.FocusOnMeshNode(2)*3);
                target = [xPos,yPos,zPos]; 
                if isnan(target(1))
                    theMeshExploded = 1
                end
                camtarget(axesList(1), target );
                axis(axesList(1), [initialCamera(1) + xPos, initialCamera(2) + xPos, initialCamera(3) + yPos, initialCamera(4) + yPos, initialCamera(5) + zPos, initialCamera(6) + zPos]);
            end

            %% -------- END OF DRAWING ------------
            % ----- Contacts -----
            ticContact = tic;
            
            cache.prevContactIDs = cache.contactIDs;
            Jc = sparse( 0, mesh3D.N*3 );
            phi = [];
            cInfos = contactInfo3D.empty;
            for j = 1:numel(contactFinders)
                [ Jci, phii, cInfosi ] = contactFinders{j}.findContacts(mesh3D, frame);
                Jc = [ Jc; Jci ];
                phi = [ phi; phii];
                cInfos = [ cInfos, cInfosi ];
            end
            cache.contactIDs = [ cInfos(:).contactID ];
            cache.cInfo = cInfos;
            td{k}.contactCount = numel(phi); 
            td{k}.lastContact = toc(ticContact);

            %------ execute the scripted animations
            mesh3D.animationDOFs = [];
            for j = 1:numel(animationScripters)
                animationScripters{j}.scriptMesh(mesh3D, integrator, frame,h);
            end

            % ----- full F C and Forces needed for Elastification -----
                    % (and subsets of these will be used in the adaptive rigid
                    % solve)
            
            ticlastSimulate = tic; % last simulate does not include the collision detection, only the time integration, elastification and rigidification
            
            startComputeFCForces = tic;
            cache.oldF = cache.F;
            
            cache.F = mesh3D.B * mesh3D.p;
            
            mesh3D.animationInds = floor(mesh3D.animationDOFs./3);

            energyModel.computeEnergy(mesh3D, cache.F)
            cache.C = energyModel.derivative2HessianC;
            cache.dpsidF = energyModel.derivative1Gradient;
            cache.elasticForces = energyModel.elasticForces; 

            td{k}.fullFCForces = toc( startComputeFCForces );

            % ----- RIGIDIFY -----
            cache.cInfo = cInfos;

            % ----- ELASTIFICATION -----

            if isa( mesh3D, 'AdaptiveMesh3D' )
                ticElastifyQuickSolve = tic;
                quickSolve3D( cache, integrator, mesh3D, h, Jc, phi, settings, animationScripters, frame); 
                quicksolveTime = toc(ticElastifyQuickSolve);
                td{k}.lastElastifyQuickSolve = quicksolveTime;
                
                ticElastify = tic;
                rigidificator.checkForElastification( mesh3D, cache, frame, h, settings );
                td{k}.lastElastification = toc(ticElastify);
            else
                td{k}.lastElastifyQuickSolve = 0;
                td{k}.lastElastification = 0;
            end

            if ( settings.PlotEDotHist )
                EDotNorms = cache.edotnorms;
                histogram( axesList(2), log10( EDotNorms ), [-inf -7:.1:1 inf], 'Normalization', 'probability' );
                title(axesList(2),'log10 EDot Fro Norm Squared');
                ylim( axesList(2), [0,0.1] );
                xline( axesList(2), log10(rigidificator.RigidificationThreshold) );

                EDotApproxNorms = cache.edotApproxNorms;
                histogram( axesList(3), log10( EDotApproxNorms ), [-inf -7:.1:1 inf], 'Normalization', 'probability' );
                title(axesList(3),'log10 EDotApprox Fro Norm Squared');
                ylim( axesList(3), [0,0.1] );
                xline( axesList(3), log10(rigidificator.ElastificationThreshold) );
            end

            % ------- INTEGRATION -------------
            % Advance the state of the meshes
            if ~settings.quicksolveSimulation
                integrator.integrate(mesh3D, h, Jc, phi, cInfos, settings, cache, td{k}, animationScripters, frame);
            end
            elapsed = elapsed + h;
            
            td{k}.lastSimulate = toc( ticlastSimulate );
            % ------ EXTRA DATA COLLECTION -----
            if settings.RecordFramePositionInTD 
                td{k}.p = [td{k}.p,mesh3D.p];
            end
            
            countParticles = 0;
            countElements = 0;    
            countTotalParticles = 0;
            countTotalTris = 0;

            countTotalParticles = countTotalParticles + mesh3D.N;
            countTotalTris = countTotalTris + size(mesh3D.t, 1);
            
            countTotalDofs = 0;
            
            if isa(mesh3D, 'AdaptiveMesh3D')
               countParticles = countParticles + numel(mesh3D.ElasticInds);
               countElements = countElements + numel(mesh3D.ElasticTetInds);
            else
               countParticles = countParticles + mesh3D.N;
               countElements = countElements + size(mesh3D.t, 1);
            end

            td{k}.countParticles = countParticles;
            td{k}.countTris = countElements;
            td{k}.countTotalParticles = countTotalParticles;
            td{k}.countTotalTris = countTotalTris;
            rigidBodies = 0;

            if isa(mesh3D, 'AdaptiveMesh3D')
                rigidBodies = rigidBodies + numel(mesh3D.RigidBodies);
            end
            td{k}.rigidBodies = rigidBodies;
            
            totalDofs = 0;
            td{k}.totalDofs = totalDofs + numel(mesh3D.activeDOFs);

            if isa(mesh3D, 'AdaptiveMesh3D')
                td{k}.linearMomentum = [ td{k}.linearMomentum, ...
                                      [ mesh3D.getAdaptiveLinearMomentum; ...
                                        mesh3D.getLinearMomentum ] ];
            else
                td{k}.linearMomentum = [ td{k}.linearMomentum, ...
                                      [ mesh3D.getLinearMomentum;
                                        mesh3D.getLinearMomentum ] ];            
            end

            td{k}.stepSum = h*frame;
            td{k}.logData;
            if settings.PlotMomentum
                frames = 1:size(td{k}.linearMomentum,2);
                plot( axesList(2), td{k}.linearMomentum(1:3,frames)' );
                title( axesList(2), 'linear momentum');
                plot( axesList(3), td{k}.logCounts(6,frames) )
                hold( axesList(3), 'on' );
                plot( axesList(3), td{k}.logCounts(5,frames)-td{k}.logCounts(3,frames) )
                plot( axesList(3), td{k}.logCounts(4,frames)-td{k}.logCounts(2,frames) )
                title( axesList(3), 'bodies, rtris, rdofs' );
                hold ( axesList(3), 'off' );
            end

            % Draw the mouse
            if ( mouseLine == 0 )
                mp = [nan,nan,nan];
                if ~isnan( mouseGrabMeshId )
                    mp =  mesh3D.p(mouseGrabMeshId*3-2:mouseGrabMeshId*3)';
                end
                mouseLine = line( [mp(1), mp(1) + pullDirection(1)], [mp(2), mp(2) + pullDirection(2)], [mp(3), mp(3) + pullDirection(3)]);
            else
                mp = [nan,nan,nan];
                if ~isnan( mouseGrabMeshId )
                    mp =  mesh3D.p(mouseGrabMeshId*3-2:mouseGrabMeshId*3)';
                end
                mouseLine.XData = [mp(1), mp(1) + pullDirection(1)];
                mouseLine.YData = [mp(2), mp(2) + pullDirection(2)];
                mouseLine.ZData = [mp(3), mp(3) + pullDirection(3)];
            end
        end
    end

    tEnd  = cputime - tStart;
    disp("cpuTime for this run: " + tEnd + " seconds" );
    td{1}.SimulationLoopFull = tEnd;

    if settings.MakeVideo == 1
        close(v);
    end

    % if a snapshot doesn't exist, save an image of what just ran
    [ST,~] = dbstack('-completenames');
    [pathstr, name, ~] = fileparts( ST(2).file );
    jpgpath = strcat( pathstr, filesep, name, ".jpg" );
    if ~isfile(jpgpath)
        if isvalid( mainFig )
            f = getframe(mainFig);
            imwrite( f.cdata, jpgpath );
        else
            disp('quit with ESC to save image for visual table of contents');
        end
    end

    %% -------------------------------------------------------

    function onKeyPressed(~, event)
        if strcmp( event.Key, 'r' )
            % see the 2D version for inspiration on how to implement
            for i2 = 1:comparisons
                meshes{i2} = initialMeshes{i2}.clone();
                caches{i2}.ActiveB = meshes{i2}.B;
                caches{i2}.clear();  % clear without deleting precomputation
            end
            
            for i2 = 1:numel(animationScripters)
                animationScripters{i2}.reset();
            end
            frame = 0;
            elapsed = 0;
            
            cfs = contactFinders; %{i3};
            for j3 = 1:numel(cfs)
                cfs{j3}.plotHandle = 0;
            end
            for k2 = 1:numel(td)
                td{k2}.log = [];
                td{k2}.logCounts = [];
                td{k2}.linearMomentum = [];
            end
            cla;
            if settings.Shading
                camlight(axesList(1),'right', 'infinite');
                shading(axesList(1),'interp'); % interp;
                lighting(axesList(1),'gouraud'); % gouraud;
            end
            mouseDown = 0;
            mousePoint = [nan, nan, nan]; % Current mouse point
            pullDirection = [nan nan nan];
            mouseGrabMeshId = nan;  % Grabbed mesh
            mouseLine = 0;
            disp('Simulation reset');
        elseif strcmp(event.Key, 'p') || strcmp(event.Key, 'space') 
            if running == 1
                running = 0;
                disp('Simulation paused');
            else 
                running = 1;
                disp('Simulation started');
            end
        elseif strcmp(event.Key, 's')
           stepOnce = 1;
        elseif strcmp(event.Key, 'd')
           settings.RigidificationEnabled = ~settings.RigidificationEnabled;
           disp( "Rigidification = " + settings.RigidificationEnabled );
        elseif strcmp(event.Key, 'l')
            settings.DrawLambdas = ~settings.DrawLambdas;
            disp( "DrawLambdas = " + settings.DrawLambdas );
        elseif strcmp(event.Key, 'e')
            settings.ElastificationEnabled = ~settings.ElastificationEnabled;
            disp( "elastification = " + settings.ElastificationEnabled );
        elseif strcmp(event.Key, 't')
            rigidificators.updateRigidBodies(meshes, [1,3,10,16,44, 29, 57, 14, 39, 70, 88, 100, 69,79,99, 160], settings);
        elseif strcmp( event.Key, 'escape') 
            % remove handlers from figure
%             set(mainFig, 'WindowButtonUpFcn', '' );
%             set(mainFig, 'WindowButtonDownFtcn', '' );
%             set(mainFig, 'WindowButtonMotionFcn', '' );
             set(mainFig, 'WindowKeyPressFcn', '' );
            notDone = false;
        end        
    end

    function onMouseDrag(~, event)
        if mouseDown && (~isnan(mouseGrabMeshId))
            mousePos = get(gca,'CurrentPoint');
            R = getCameraTransform(gca);
            mousePoint = R*[mousePos(1,:)]'; %projects the mouse point
            mousePoint = mousePoint(1:3)';

            pointDOFs = mouseGrabMeshId*3-2:mouseGrabMeshId*3;
            vertexPointPos = meshes{1}.p(pointDOFs)';
            pullDirection = (vertexPointPos-mousePoint)';
            for i2 = 1:comparisons
                constMultiplyer = 1;
                if ~meshes{1}.pinned(mouseGrabMeshId)
                    meshes{i2}.v(pointDOFs) = meshes{i2}.v(pointDOFs) + constMultiplyer * pullDirection;
                end
            end
        end
    end

    function R = getCameraTransform(figure)
        camPosition = get(figure, 'CameraPosition'); % camera position
        camTarget = get(figure, 'CameraTarget'); % where the camera is pointing to

        camDirection = camPosition - camTarget; % camera direction
        camUp = get(figure, 'CameraUpVector'); % camera 'up' vector

        % build an orthonormal frame based on the viewing direction and the 
        % up vector (the "view frame")
        zAxis = camDirection/norm(camDirection);    
        upAxis = camUp/norm(camUp); 
        xAxis = cross(upAxis, zAxis);
        yAxis = cross(zAxis, xAxis);

        R = [xAxis; yAxis; zAxis];
    end

    function onMouseDown(~, event)
        %TODO identify triangle and coords
        % start grabbing and on mouse move apply force
        % mouse release stop grabbing
        
        % Note seems funny that the position is not available from the
        % argumetns, but this seems to be the way matlab does it!
        R = getCameraTransform(gca);
        mousePos = get(gca,'CurrentPoint');
        mousePoint = R*[mousePos(1,:)]'; %projects the mouse point
        mousePoint = mousePoint(1:3)'; %only take the first 3 coordinates, no need for w.
        plot3(mousePoint(1),mousePoint(2),mousePoint(3),'.');

        mouseGrabMeshId = NaN;
        %Don't need to iterate comparisons as we use the same tet
        for i2 = 1:comparisons
            mesh3D = meshes{i2};
            if ( ~isnan( mouseGrabMeshId ) )
                continue; % already found something... 
            end

            P = mesh3D.getPositionFormatted;
            nearestPointInd = dsearchn(P,mousePoint);
            mouseGrabMeshId = nearestPointInd;
            
        end  
        mouseDown = 1;
    end

    function onMouseUp(~, event)
        mouseDown = 0;
        mouseGrabMeshId = nan;
    end

end