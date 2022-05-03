function [td] = simulate( meshes, integrators, h, settings, rigidificators, contactFinders, meshContactFinder, animationScripter)
    %SIMULATE Runs a single or multiple simulations and displays them into
    %a figure
    %   Each parameter can either be specified as itself or as a cell array
    %   to support both multi simulations (simulation comparisons) or
    %   single simulations.
    %
    %   You can interact with the simulation while its running by pressing
    %   p to play/pause, r to reset and s to step forward in time. Clicking
    %   anywhere in the simulation will give an impulse to the closest
    %   particle.
    %
    %
    % Inputs:
    %   meshes: a cell array of (or a single) mesh that specify the 
    %       simulation mesh.  If you want multiple meshes, then merge them!
    %   integrators: a cell array of (or a single) integrator (see Integrator 
    %       class) that indicates how to step forward the simulation
    %   h: time step of the simulation
    %   settings: settings of the simulation, see class SimulationSettings
    %   rigidificators: a cell array of (or a single) rigidificator that
    %       implements the rigidification and unrigidification process
    %   contactFinders: a cell array of (or a single) contact finder that
    %       implements the process in which to find contacts and returns
    %       constraints for those contacts
    %   mcf: the mesh self collision contact detection object
    
    % Note cyan variable colouring indicates variables used within several
    % functions defined within this function, i.e., the scope spans
    % multiple functions!
    
    if nargin < 8
        animationScripter = {NullAnimationScript()};
    end
    
    if nargin < 7
        %no inter elastic contact finder by default
        meshContactFinder = NullContactFinder();
    end    
    if nargin < 6
        contactFinders = NullContactFinder();
    end
    if nargin < 5
        rigidificators = EDotMexRigidificator();
    end
    if nargin < 4
       settings = SimulationSettings(); 
    end
    if nargin < 3
        h = 0.01; 
    end
    if nargin < 2
        integrators = LDLBackwardEuler();
    end
    if nargin < 1
        warning( 'MUST provide at least the mesh to simulate' );
        return;
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
    
    if ~isa(animationScripter, 'cell')
        animationScripter = { animationScripter };
    end
    
   [comparisons, initialMeshes, meshes, integrators, rigidificators] = setupComparisons2D(meshes,integrators,rigidificators);

    [mainFig,axesList, initialCamera] = setupWindow2D(settings, meshes, integrators);

    caches = setupCache2D(settings, meshes, rigidificators, comparisons, h);

    % prep for video generation 
    video = 0;
    [ST,~] = dbstack('-completenames');
    [~, exampleName, ~] = fileparts( ST(2).file );
    if settings.MakeVideo
        stamp = string(datetime('now','Format','y-MMM_d-HH_mm_ss'));
        vidFileName = strcat( 'out', filesep, exampleName, stamp, ".mp4" );
        video = VideoWriter(vidFileName,'MPEG-4');
        open(video);
    end

    %events
    set(mainFig, 'Name', 'Simulation');
    set(mainFig, 'WindowKeyPressFcn', @onKeyPressed);
    set(mainFig, 'WindowButtonUpFcn', @onMouseUp);
    set(mainFig, 'WindowButtonDownFcn', @onMouseDown);
    set(mainFig, 'WindowButtonMotionFcn',@onMouseDrag);
    set(mainFig, 'color', 'w');

    plots = cell( comparisons, 1 );
    names = cell( comparisons, 1 );
    legendSet = 0;
    % need init first??
    simulTimeText.String        = "Simulation time:";
    realTimeText.String         = "Wall clock time:";
    avgFpsText.String           = "Avg FPS:"; 
    renderTimeText.String       = "Render time:";
    elastificationText.String   = "Elastification time: ";
    rigidificationText.String   = "Rigidification time: ";
    contactText.String          = "Contact detect time: ";
    simulateText.String         = "Simulate time:       ";
    particleText.String         = "Elastic particles: ";
    triangleText.String         = "Elastic triangles: ";
    rigidBodyText.String        = "Rigid bodies: ";
    
    printKeyboardControls;
        
    %UI variables
    mouseDown = 0;
    closestTrianglePositionID = 0;
    closestMesh = NaN;
    mousePoint = [NaN NaN]; % Current mouse point
    mouseGrabMeshId = NaN;  % Grabbed mesh
    mouseGrabTriInd = NaN;  % Grabbed triangle
    mouseGrabBCC = [1,0,0]; % Barycentric coordaintes of the grabbed point
    mouseLine = 0;
    
    %timing
    elapsed = 0;
    realElapsed = 0;
    frame = 0;
    running = settings.RunRightAway;
    stepOnce = 0;
    lastRender = 0;
    td = cell( comparisons, 1 );
    for k = 1:comparisons
        td{k} = TimingData();
    end

    singleUseFigures = [];
    
    Display();
    % initialize the rendering (text, legend, etc.)
    init();
    drawingNeedsUpdating = false;
    
    notDone = true;
    
    skip = settings.PlotSkip + 1;
    while ishghandle(mainFig) && notDone

        % If we are not running, or don't have a step request, then we
        % really shouldn't bother doing all the redrawing!  Only collect
        % video frames for when the simulation is running!
        if ~running
            if stepOnce
                disp('Stepping once');
                stepOnce = 0; 
            elseif drawingNeedsUpdating
                DrawAndDisplay;
                continue;
            else
                pause(0.2);
                continue;
            end
        end 
        
        frame = frame + 1;
        
        ticFrame = tic;

        if frame >= settings.FramesToRecord && settings.FramesToRecord ~= -1
            notDone = false;
        end
        
        %% -------- BEGIN SIMULATION CODE -----
       
        for k = 1:comparisons
            mesh2D = meshes{k};
            rigidificator = rigidificators{k};

            contactFinderList = contactFinders; % No longer part of the comparisons
            integrator = integrators{k};
            cache = caches{k};
            td{k}.lastSimulate = tic;
            
            % ----- COLLISION DETECTION -----
            ticContact = tic;
            
            cache.prevContactIDs = cache.contactIDs;
            
            Jc = []; % need to set this to the right size? might just work like this??
            phi = [];
            cInfos = contactInfo.empty;
            for i = 1:numel(contactFinderList)
                [ Jci, phii, cInfosi ] = contactFinderList{i}.findContacts( mesh2D, elapsed );
                Jc = [ Jc; Jci ];
                phi = [ phi; phii];
                cInfos = [ cInfos, cInfosi ];
            end
            %find elastic contacts and add them to contacts
            if settings.ElasticContacts
                [ JcElastic, phiElastic, cInfosOut ] = meshContactFinder.findContacts( mesh2D );
                Jc = [ Jc; JcElastic ];
                phi = [ phi; phiElastic ];
                cInfos = [ cInfos, cInfosOut ];
            end
            cache.contactIDs = [ cInfos(:).contactID ];
            
            td{k}.contactCount = numel(phi); 
            td{k}.lastContact = toc(ticContact);
            
            % ----- full F C and Forces needed for Elastification -----
            % (and subsets of these will be used in the adaptive rigid
            % solve)
            
            startComputeFCForces = tic;
            
            %execute the scripted animations
            mesh2D.animationDOFs = [];
            for j = 1:numel(animationScripter)
                animationScripter{j}.scriptMesh(mesh2D, integrator, frame, h);
            end
            mesh2D.animationInds = floor((mesh2D.animationDOFs+1)./2);
            
            cache.oldF = cache.F;
            
            cache.F = mesh2D.B * mesh2D.p;
            
            [ii, jj, CblockVals, cache.dpsidF] = mexComputeSTVKGradHess2D( cache.F, mesh2D.elA, mesh2D.elMu, mesh2D.elLambda );
            cache.C = sparse( ii, jj, CblockVals );
            cache.elasticForces = mesh2D.B' * cache.dpsidF; 
            cache.elasticForces(mesh2D.pinnedDOFs) = 0;

            td{k}.fullFCForces = toc( startComputeFCForces );
            
            ticSimulate = tic;

            % ----- RIGIDIFY ----- 
            cache.cInfo = cInfos;
            isAdaptiveComparision = isa( mesh2D, 'AdaptiveMesh' );

            if isAdaptiveComparision
                ticElastifyQuickSolve = tic;
                cache.cInfo = cInfos;
                quickSolveNoConstraint( cache, integrator, mesh2D, h, Jc, phi, settings, animationScripter, frame); 
                td{k}.lastElastifyQuickSolve = toc(ticElastifyQuickSolve);
                
                ticElastify = tic;
                rigidificator.checkForElastification( mesh2D, cache, frame, h, settings );
                td{k}.lastElastification = toc(ticElastify);
            else
                td{k}.lastElastifyQuickSolve = 0;
                td{k}.lastElastification = 0;
            end

            if ( settings.PlotEDotHist && isa(mesh2D,"AdaptiveMesh")&& (k == 1 || k == 2))
                EDotNorms = cache.edotnorms;
                histogram( axesList(2), log10( EDotNorms ), [-inf -10:.1:1 inf], 'Normalization', 'probability' );
                title(axesList(2),'log10 EDot Fro Norm Squared');
                ylim( axesList(2), [0,0.1] );
                xline( axesList(2), log10(rigidificator.RigidificationThreshold) );   
                EDotApproxNorms = cache.edotapproxnorms;
                histogram( axesList(3), log10( EDotApproxNorms ), [-inf -7:.1:1 inf], 'Normalization', 'probability' );
                title(axesList(3),'log10 EDotApprox Fro Norm Squared');
                ylim( axesList(3), [0,0.1] );
                xline( axesList(3), log10(rigidificator.ElastificationThreshold));
            elseif settings.PlotPhiHist
                plot( axesList(2), log(-phi) );
                 title(axesList(2),'log phi of contacts');
                 ylim( axesList(2), [-8,-2] );
            elseif settings.PlotApproxRateWork
                rateOfWork = zeros(size(mesh2D, 1),1);
                vel = mesh2D.v + cache.ApproximatedDeltaV; 
                F =  mesh2D.B * (mesh2D.p + h*vel+mesh2D.p)/2;
                FDot =  mesh2D.B * vel;
                for ii = 1:size(mesh2D.t, 1)
                    stress = -cache.dpsidF(4*ii-3:4*ii)./mesh2D.el(ii).area;
                    stress = reshape(stress,2,2);
                    localFDot = reshape(FDot(4*ii-3:4*ii),2,2);
                    localF = reshape(F(4*ii-3:4*ii),2,2);
                    strainRate = 0.5 * (localFDot' * localF + localF' * localFDot);
                    workSquared = (stress.*strainRate).^2;
                    rateOfWork(ii) = sum(workSquared(:));
                end
                plot( axesList(2), log(rateOfWork));
                title(axesList(2),'log approx rate of work');
                ylim( axesList(2), [-50,50] );
            elseif settings.PlotPolarDecomposition
                a = gcf;
                for i = 1:1%mesh.N
                    smallF = reshape(cache.F(4*i-3:4*i),2,2);
                    [R,U,V] = polarDecomposition(smallF);

                    hold on;
                    imagesc(axesList(2), U);
                    title(axesList(2),'Unitary');
                    hold on;
                    imagesc(axesList(3), R);
                    title(axesList(3),'Hermitian');
                end
            end
            
            
            % ------- INTEGRATION -------------
            if ~settings.quicksolveSimulation
                integrator.integrate( mesh2D, h, Jc, phi, cInfos, settings, cache, td{k}, animationScripter, frame ); 
            end
            % ----------------------
            td{k}.lastSimulate = toc( ticSimulate );
            td{k}.logData;
            
            if settings.RecordFramePositionInTD 
                td{k}.p = [td{k}.p,mesh2D.p];
            end
            
            %%%% ---- print timings -----
            if ( settings.PrintTimings )
                disp("Timings for comparision " + k );
                td{k}
            end
            
        end
        
        elapsed = elapsed + h;
        realH = toc(ticFrame);
        realElapsed = realElapsed + realH;
        
        % Collect other timing data related quantities
        for k = 1:comparisons
            countParticles = 0;
            countTris = 0;    
            countTotalParticles = 0;
            countTotalTris = 0;
            totalDofs = 0;
            
            countTotalParticles = countTotalParticles + meshes{k}.N;
            countTotalTris = countTotalTris + size(meshes{k}.t, 1);
            if isa(meshes{k}, 'AdaptiveMesh')
               countParticles = countParticles + numel(meshes{k}.ElasticInds);
               countTris = countTris + numel(meshes{k}.ElasticTriInds);
               totalDofs = totalDofs + numel(meshes{k}.activeDOFs);
            else
               countParticles = countParticles + meshes{k}.N;
               countTris = countTris + size(meshes{k}.t, 1);
               totalDofs = totalDofs + numel(meshes{k}.activeDOFs);
            end

            td{k}.countParticles = countParticles;
            td{k}.countTris = countTris;
            td{k}.countTotalParticles = countTotalParticles;
            td{k}.countTotalTris = countTotalTris;
            td{k}.totalDofs = totalDofs;
            rigidBodies = 0;

            mesh2D = meshes{k};
            if isa(mesh2D, 'AdaptiveMesh')
                rigidBodies = rigidBodies + numel(mesh2D.RigidBodies);
            end
            td{k}.rigidBodies = rigidBodies;
                       
            mesh2D = meshes{k};
            if isa(mesh2D, 'AdaptiveMesh')
                td{k}.linearMomentum = [ td{k}.linearMomentum, ...
                    [ mesh2D.getAdaptiveLinearMomentum; ...
                      mesh2D.getLinearMomentum ] ];
            end
             
        end
         
        %% -------- END SIMULATION CODE -----
        %% -------- BEGIN DRAWING -----
        DrawAndDisplay;

        if ( settings.PlotContactForceProfile )
            points = [cache.cInfo.point];
            xp = points(1:2:end);
            if numel(xp) == 0
                cla(axesList(2))
            else
                hold(axesList(2),'off'); 
                plot(axesList(2),xp,cache.prevLambdas(1:2:end),'r')
                hold(axesList(2),'on'); 
                plot(axesList(2),xp, cache.prevLambdas(2:2:end),'b');
                %legend(axesList(2),'normal force', 'friction force');
                % legend doesn't work not sure why
                xlabel(axesList(2),'contacts x position');
                ylabel(axesList(2),'force');
            end
        elseif settings.PlotRateWork
                rateOfWork = zeros(size(mesh2D, 1),1);
                F =  mesh2D.B * mesh2D.p;
                FDot =  mesh2D.B * mesh2D.v;
                for ii = 1:size(mesh2D.t, 1)
                    stress = -cache.dpsidF(4*ii-3:4*ii)./mesh2D.el(ii).area;
                    stress = reshape(stress,2,2);
                    localFDot = reshape(FDot(4*ii-3:4*ii),2,2);
                    localF = reshape(F(4*ii-3:4*ii),2,2);
                    strainRate = 0.5 * (localFDot' * localF + localF' * localFDot);
                    workSquared = (stress.*strainRate).^2;
                    rateOfWork(ii) = sum(workSquared(:));
                end
                plot( axesList(2), log(rateOfWork));
                title(axesList(2),'log rate of work');
                ylim( axesList(2), [-50,50] );
        end
        %% -------- END OF DRAWING ------------
    end
    
    disp('Simulation finished');

    if settings.MakeVideo
        close(video);
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

    function DrawAndDisplay()
        drawingNeedsUpdating = false;
        
        delete(singleUseFigures);
        singleUseFigures = [];

        if mod(frame, skip) == 0
            Display();

            % ---- Draw collision objects if necessar (i.e., fi mobile)
            %for k = 1:comparisons
                cfs = contactFinders;%{k};
                for j = 1:numel(cfs)
                    cfs{j}.render( elapsed );
                end
            %end

            % -------- DRAW CONTACTS ---------
            for k = 1:comparisons
                cache = caches{k};
                cInfos = cache.cInfo;

                if (settings.DrawContact || settings.DrawLambdas)
                    if ( isempty(cInfos) ) 
                        continue;
                    end                
                    for j = 1:numel(cInfos)
                        cp = cInfos(j).point;
                        if settings.DrawContact
                            f = scatter( cp(1), cp(2), 35, 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5 );
                            singleUseFigures = [singleUseFigures,f];
                            if ~isnan(cInfos(j).pointAlpha)
                                cp2 = cInfos(j).pointAlpha;
                                f = scatter( cp2(1), cp2(2), 35, 'MarkerFaceColor', [0, 1, 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
                                singleUseFigures = [singleUseFigures,f];
                                
                                cp2 = cInfos(j).pointAlphaInv;
                                f = scatter( cp2(1), cp2(2), 35, 'MarkerFaceColor', [1, 0, 1], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
                                singleUseFigures = [singleUseFigures,f];
                            end
                        end
                        if settings.DrawLambdas
                            shorten = settings.DrawLambdasScale;
                            cn = cInfos(j).normal;
                            ct = cInfos(j).tangent;
                            cl = cache.prevLambdas(j*2-1:j*2);
                            cf = (cn*cl(1) + ct*cl(2)) * shorten;
                            f = line( [cp(1), cp(1)-cf(1)], [cp(2), cp(2)-cf(2)]);
                            singleUseFigures = [singleUseFigures,f];
                        end
                    end
                end
            end

            % Draw the mouse
            
            if ( mouseLine == 0 )
                mp = [nan,nan];
                if ~isnan( mouseGrabMeshId )
                    mpts = reshape(meshes{1}(mouseGrabMeshId).p,2,meshes{1}(mouseGrabMeshId).N)';
                    mp = mouseGrabBCC * mpts( meshes{1}(mouseGrabMeshId).t(mouseGrabTriInd), : );
                end
                mouseLine = line( [mp(1), mousePoint(1)], [mp(2), mousePoint(2)]);
            else
                mp = [nan,nan];
                if ~isnan( mouseGrabMeshId )
                    mpts = reshape(meshes{1}(mouseGrabMeshId).p,2,meshes{1}(mouseGrabMeshId).N)';
                    triPts = mpts( meshes{1}(mouseGrabMeshId).t(mouseGrabTriInd,:), : );
                    mp = mouseGrabBCC * triPts;
                end
                mouseLine.XData = [mp(1), mousePoint(1)];
                mouseLine.YData = [mp(2), mousePoint(2)];
            end


            if settings.MakeVideo
                F = getframe(mainFig);
                writeVideo(video, F);
            else
                drawnow; % seems like getframe does it
            end
        end
    end
   
    function Display()
        
        ticRender = tic;
        set(0, 'CurrentFigure', mainFig)
        
        % retains plots in the current axes so that new plots added to the axes do not delete existing plots
        hold on;
        
        if ~isempty(settings.FocusOnMeshNode)
            meshesFocus = meshes{settings.FocusOnMeshNode(1)};
            xPos = meshesFocus.p(settings.FocusOnMeshNode(2)*2-1);
            yPos = meshesFocus.p(settings.FocusOnMeshNode(2)*2);
            heightadjustment = (initialCamera(4)-initialCamera(3))/2.0;
            axis(gca, [initialCamera(1) + xPos, initialCamera(2) + xPos, initialCamera(3) + yPos + heightadjustment, initialCamera(4) + yPos + heightadjustment]);
        end
        
        % plot the meshes for each comparisons
        for kk = 1:comparisons
            singleUseFigures = plotMesh( h, meshes{kk}, settings, singleUseFigures, caches{kk}, integrators{kk}, rigidificators{kk});

            % alpha = 1 if a single simulation
            % only override the alpha if we are doing comparisons?
            if ( comparisons > 1 ) 
                %alpha(1 / (comparisons * comparisons));
            end
        end   
        % update the legend if needed
        if settings.DrawLegend && ~legendSet && comparisons > 1
            for kk = 1:comparisons
                plots(kk) = bar(0, 0, 'FaceColor', genColor(k));
                names{kk} = integrators{k}.Name;
            end
            legend(plots, names);
            legendSet = 1;
        end 
        
        % update the text
        if settings.DrawTimings
            simulTimeText.String = sprintf('Simulation time: %f', elapsed);
            realTimeText.String  = sprintf('Wall clock time: %f', realElapsed);
            avgFpsText.String = sprintf('Avg FPS: %i', floor(frame / realElapsed));
            renderTimeText.String       = sprintf('Render time:         %f', lastRender);
            elastificationText.String = "Elastification time: " + ...
                join(cellfun(@(x) x.lastElastification + " (" + x.lastElastifyQuickSolve + " quick solve)", td ));
            rigidificationText.String   = "Rigidification time: " + join(cellfun(@(x) ""+x.lastRigidification, td ));
            contactText.String          = "Contact detect time: " + ...
                join(cellfun(@(x) x.lastContact +" ("+ x.contactCount+")", td ));
            simulateText.String         = "Simulate time:       " + join(cellfun( @(x) ""+x.lastSimulate, td ));
            particleText.String = "Elastic particles: " + join(cellfun(@(x) x.countParticles + " (out of " + x.countTotalParticles + ")", td ));
            triangleText.String = "Elastic triangles: " + join(cellfun(@(x) x.countTris + " (out of " + x.countTotalTris + ")", td ));
            rigidBodyText.String = "Rigid bodies: " + join(cellfun(@(x) ""+x.rigidBodies, td ));
        end
        lastRender = toc(ticRender);
    end

    function printKeyboardControls
        disp("-----------------------------------------");
        disp("Keyboard Controls (type in figure window)  ");
        disp("-----------------------------------------");
        disp("Basic Controls");
        disp("    r            reset                         ");
        disp("    p space      start and stop simulation     ");
        disp("    s            step once                     ");
        disp("    escape       quit simulation               ");
        disp("    h            help! this menu               ");
        disp("Visulaization");
        disp("    c            toggle draw contact locations ");
        disp("    l            toggle draw contact forces");
        disp("    , .          smaller/bigger contact force scale");        
        disp("    f            toggle draw elastic forces ");
        disp("    v            toggle draw velocities   ");
        disp("    k            toggle draw approx velocities   ");
        disp("    m            toggle draw approx EDots   ");
        disp("Rigidification");
        disp("    e            toggle elastification    ");
        disp("    d            toggle rigidification    ");
        disp("    u            unrigidify all           ");
        disp("    z            zero out velocities      ");
        disp("    ; '          smaller/bigger tau_E     ");
        disp("    [ ]          smaller/bigger tau_R     ");
        disp("    = -          decrease/increase quicksolve PCG iterations");
        disp("Experimental");
        disp("    t            test rigidification (see code for element set) ");
        disp("    j            toggle regularizationWarmstart ?");
        disp("    1-9          rigidify numbered element ?");
        disp("-----------------------------------------");
    end

    function onKeyPressed(~, event)
        
        % not genearlly true that we need to update, but in may cases we
        % do, so easy to just always request a redraw on keypress,
        % regardless.
        drawingNeedsUpdating = true;


        if strcmp(event.Key, 'tab')
            if ( settings.MakeVideo == 1 ) 
                disp('video recording finished');
                close(video);
                settings.MakeVideo = 0;
            else
                stamp = string(datetime('now','Format','y-MMM-d_HH-mm-ss'));
                vidFileName = strcat( 'out', filesep, exampleName, stamp, '.mp4' );
                disp(strcat("video recording started, writing to ", vidFileName));
                video = VideoWriter(vidFileName,'MPEG-4');
                open(video);
                settings.MakeVideo = 1;
            end
        elseif strcmp(event.Key, 'period')
            settings.DrawLambdasScale = settings.DrawLambdasScale*2;
        elseif strcmp(event.Key, 'comma')
            settings.DrawLambdasScale = settings.DrawLambdasScale/2;
        elseif strcmp(event.Character, '[')
            for i2 = 1:comparisons
                rigidificators{i2}.RigidificationThreshold = rigidificators{i2}.RigidificationThreshold / 10^(1/4);
                disp( "comparison " + i2 + " tau_R = " + rigidificators{i2}.RigidificationThreshold );
            end
        elseif strcmp(event.Character, ']')
            for i2 = 1:comparisons
                rigidificators{i2}.RigidificationThreshold = rigidificators{i2}.RigidificationThreshold * 10^(1/4);
                disp( "comparison " + i2 + " tau_R = " + rigidificators{i2}.RigidificationThreshold );
            end
        elseif strcmp(event.Key, 'semicolon')
            for i2 = 1:comparisons
                % constant chosen such that 3 steps is one power of 10
                rigidificators{i2}.ElastificationThreshold = rigidificators{i2}.ElastificationThreshold / 10^(1/4);
                disp( "comparison " + i2 + " tau_E = " + rigidificators{i2}.ElastificationThreshold );
            end
        elseif strcmp(event.Key, 'quote')
            for i2 = 1:comparisons
                rigidificators{i2}.ElastificationThreshold = rigidificators{i2}.ElastificationThreshold * 10^(1/4);
                disp( "comparison " + i2 + " tau_E = " + rigidificators{i2}.ElastificationThreshold );
            end
        elseif strcmp(event.Character, '-')
            settings.PCGiterations = settings.PCGiterations - 1;
            if (settings.PCGiterations < 1)
                settings.PCGiterations = 1;
            end
            disp( "quick solve PCG iterations = " + settings.PCGiterations );
        elseif strcmp(event.Character, '=')
            settings.PCGiterations = settings.PCGiterations + 1;
            disp( "quick solve PCG iterations = " + settings.PCGiterations );
        elseif strcmp(event.Key, 'h')
            printKeyboardControls;
        elseif strcmp(event.Key, 'z')
            for i2 = 1:comparisons
                meshes{i2}.v = zeros(size(meshes{i2}.v));
            end 
        elseif strcmp(event.Key, 'c')
            settings.DrawContact = ~settings.DrawContact;
        elseif strcmp(event.Key, 'd')
           settings.RigidificationEnabled = ~settings.RigidificationEnabled;
           disp( "Rigidification = " + settings.RigidificationEnabled );
        elseif strcmp(event.Key, 'e')
            settings.ElastificationEnabled = ~settings.ElastificationEnabled;
            disp( "elastification = " + settings.ElastificationEnabled );
        elseif strcmp(event.Key, 'f')
            settings.DrawForces = ~settings.DrawForces;
        elseif strcmp(event.Key, 'l')
            settings.DrawLambdas = ~settings.DrawLambdas;
        elseif strcmp(event.Key, 'j')
            settings.regularizationWarmstart = ~settings.regularizationWarmstart;
            disp( "regularization = " + settings.regularizationWarmstart );
        elseif strcmp(event.Key, 'k')
            settings.DrawApproxDv = ~settings.DrawApproxDv;
            disp( "DrawApproxDv = " + settings.DrawApproxDv );
        elseif strcmp(event.Key, 'm')
            settings.DrawApproxEDots = ~settings.DrawApproxEDots;
            disp( "DrawApproxEDots = " + settings.DrawApproxEDots );
        elseif strcmp(event.Key, 'p') || strcmp(event.Key, 'space') 
            if running == 1
                running = 0;
                disp('Simulation paused');
            else 
                running = 1;
                disp('Simulation started');
            end
        elseif strcmp( event.Key, 'r' )
            for i2 = 1:comparisons
                meshes{i2} = initialMeshes{i2}.clone();
                caches{i2}.ActiveB = meshes{i2}.B;
                caches{i2}.clear();  % clear without deleting precomputation
            end
            
            for i2 = 1:numel(animationScripter)
                animationScripter{i2}.reset();
            end
            
            frame = 0;
            elapsed = 0;
            realElapsed = 0;
            legendSet = 0;
            mouseLine = 0;
%            for i3 = 1:comparisons
                % Need to clear the plots of the contact finders because
                % the whole figure is going to be cleared and needs to be
                % redrawn
                cfs = contactFinders; %{i3};
                for j3 = 1:numel(cfs)
                    cfs{j3}.plotHandle = 0;
                end
%            end
            cla;
            init();
            Display();
            disp('Simulation reset');
        elseif strcmp(event.Key, 's')
           stepOnce = 1;
        elseif strcmp(event.Key, 't')
%             rigidificators{1}.updateRigidBodies(meshes{1}, [1:11,43:47,21,52,53,55, 23,13], settings);
            disp( "rigidifying test set" );
            
            % leave elastic everything with y < 2.75 and x > 0
            bc = barycenter(reshape( meshes{1}.p, 2, [] )',meshes{1}.t);
            ids = find((bc(:,1) < 0 | bc(:,2) > -2.7));
            rigidificators{1}.updateRigidBodies(meshes{1}, ids, settings);
        elseif strcmp(event.Key, 'u')
            for i2 = 1:comparisons
                if isa(meshes{i2}, 'AdaptiveMesh')
                   meshes{i2}.RigidBodies = RigidBody.empty;
                   meshes{i2}.updateRigidState();
                end
            end
            disp('All meshes unrigidified'); 
        elseif strcmp(event.Key, 'v')
            settings.DrawVelocities = ~settings.DrawVelocities;
        elseif strcmp( event.Key, 'escape') 
            % remove handlers from figure
            set(mainFig, 'WindowButtonUpFcn', '' );
            set(mainFig, 'WindowButtonDownFcn', '' );
            set(mainFig, 'WindowButtonMotionFcn', '' );
            set(mainFig, 'WindowKeyPressFcn', '' );
            notDone = false;
        end

        if ( numel(event.Key) == 1 )
            if ((event.Key > '0') && (event.Key <= '9')) 
                tri = event.Key - '0';
                %tri = 1:tri;
                rigidificators{1}.updateRigidBodies(meshes{1}, tri, settings);
                disp( "rigidifying tri " + tri );
            end
        end
    end

    function onMouseDrag(~, ~)
        if mouseDown && (~isnan(mouseGrabMeshId))
            mousePos = get(gca, 'CurrentPoint');
            mousePoint = mousePos(1, 1:2);
            tri = meshes{1}.t(mouseGrabTriInd,:);
            p1ID = tri(1)*2-1:tri(1)*2;
            p2ID = tri(2)*2-1:tri(2)*2;
            p3ID = tri(3)*2-1:tri(3)*2;
            p1 = meshes{1}.p(p1ID);
            p2 = meshes{1}.p(p2ID);
            p3 = meshes{1}.p(p3ID);
            Vmouse = [p1';p2';p3'];
            bary = Vmouse .* mouseGrabBCC';
            vec = sum(bary)' - mousePoint';
            % apply this to all comparions
            for i2 = 1:comparisons
                constMultiplyer = -5;
                if ~meshes{1}.pinned(tri(1))
                    meshes{i2}.v(p1ID) = meshes{i2}.v(p1ID) + constMultiplyer * vec * mouseGrabBCC(1);
                end
                if ~meshes{1}.pinned(tri(2))
                    meshes{i2}.v(p2ID) = meshes{i2}.v(p2ID) + constMultiplyer * vec * mouseGrabBCC(2);
                end
                if ~meshes{1}.pinned(tri(3))
                    meshes{i2}.v(p3ID) = meshes{i2}.v(p3ID) + constMultiplyer * vec * mouseGrabBCC(3);
                end
            end
        end
    end

    function onMouseDown(~, ~)
        
        % start grabbing and on mouse move apply force
        % mouse release stop grabbing
        
        % Note seems funny that the position is not available from the
        % argumetns, but this seems to be the way matlab does it!
        mousePos = get(gca, 'CurrentPoint');
        mousePoint = mousePos(1, 1:2);

        closestElasticPoint = 0;
        closestDist = Inf;

        mouseGrabMeshId = NaN;
        %Don't need to iterate comparisons as we use the same tri
        for i2 = 1:comparisons
            mesh2D = meshes{i2};
            if ( ~isnan( mouseGrabMeshId ) )
                continue; % already found something... 
            end
                
            T = mesh2D.t;
            P = reshape(mesh2D.p,2,mesh2D.N)';
            TR = triangulation(T,P);
            mouseGrabTriInd = pointLocation( TR, mousePoint );
            if ( isnan( mouseGrabTriInd ) ) 
                continue;
            end
            mouseGrabMeshId = 1; 
            mouseGrabBCC = cartesianToBarycentric( TR, mouseGrabTriInd, mousePoint );
        end  
        mouseDown = 1;
    end

    function onMouseUp(~, ~)
        mouseDown = 0;
        
        mouseGrabMeshId = nan;
        mouseGrabTriInd = nan;
        mouseGrabBCC = [0,0,0];
    end

    function init()
        %INIT initializes the rendering of the figure for simulation
        
        %for i3 = 1:comparisons
            cfs = contactFinders;%{i3};
            for j3 = 1:numel(cfs)
                cfs{j3}.render( 0 );
            end
        %end
        
        for itr = 1:numel(caches)
            caches{itr}.ApproximatedDeltaV = zeros(size(meshes{itr}.p));
        end
        
        if settings.DrawTimings
            % NOTE: this text is placed with units set to points, and goes
            % from bottom up.  Not sure how to get this to go top down 
            % cleanly, but best not to fight too hard on this kind of
            % thing!

            font = 'FixedWidth';
            fontsize = 7;
            yoff = fontsize + 1;
            textypos = yoff;
            xpos = fontsize;
            simulTimeText = text(xpos, textypos, '-', 'fontname', font, 'fontsize', fontsize, 'Units', 'points', 'Interpreter', 'none');
            textypos = textypos + yoff;
            realTimeText = text(xpos, textypos, '-', 'fontname', font, 'fontsize', fontsize, 'Units', 'points', 'Interpreter', 'none');
            textypos = textypos + yoff;
            avgFpsText = text(xpos, textypos, '-', 'fontname', font, 'fontsize', fontsize, 'Units', 'points', 'Interpreter', 'none');
            textypos = textypos + yoff;
            renderTimeText = text(xpos, textypos, '-', 'fontname', font, 'fontsize', fontsize, 'Units', 'points', 'Interpreter', 'none');
            textypos = textypos + yoff*1.5;
            elastificationText = text(xpos, textypos, '-', 'fontname', font, 'fontsize', fontsize, 'Units', 'points', 'Interpreter', 'none');
            textypos = textypos + yoff;
            rigidificationText = text(xpos, textypos, '-', 'fontname', font, 'fontsize', fontsize, 'Units', 'points', 'Interpreter', 'none');
            textypos = textypos + yoff;
            contactText = text(xpos, textypos, '-', 'fontname', font, 'fontsize', fontsize, 'Units', 'points', 'Interpreter', 'none');
            textypos = textypos + yoff;
            simulateText = text(xpos, textypos, '-', 'fontname', font, 'fontsize', fontsize, 'Units', 'points', 'Interpreter', 'none');
            textypos = textypos + yoff*1.5;

            particleText = text(xpos, textypos, '-', 'fontname', font, 'fontsize', fontsize, 'Units', 'points', 'Interpreter', 'none');
            textypos = textypos + yoff;

            triangleText = text(xpos, textypos, '-', 'fontname', font, 'fontsize', fontsize, 'Units', 'points', 'Interpreter', 'none');
            textypos = textypos + yoff;

            if isa(meshes, 'AdaptiveMesh')
                rigidBodyText = text(xpos, textypos, '-', 'fontname', font, 'fontsize', fontsize, 'Units', 'points', 'Interpreter', 'none');
            end 
        end

        if settings.DrawLegend
            plots = zeros(1, comparisons);
            names = cell(1, comparisons);
            legendSet = 0;
        end
    end
end


