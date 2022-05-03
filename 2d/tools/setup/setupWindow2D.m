function [mainFig,axesList, initialCamera] = setupWindow2D(settings,meshes, integrators)
% [mainFig,axesList, initialCamera] = setupWindow2D(settings,meshes, integrators)
% sets up the figure window to plot the scene with respect to the settings.
    minX = 0;
    minY = 0;
    maxX = 0;
    maxY = 0;

    totalMeshCount = 0;
    
    % prepare the meshes, clone them for reset and find camera bounds
    for k = 1:numel(meshes)
        minX = min(minX, min(meshes{k}.p(1:2:end)));
        minY = min(minY, min(meshes{k}.p(2:2:end)));
        maxX = max(maxX, max(meshes{k}.p(1:2:end)));
        maxY = max(maxY, max(meshes{k}.p(2:2:end)));
        meshes{k}.prepare(integrators{k}.InfMassForPinned);
        totalMeshCount = totalMeshCount + 1;
    end
    
    % apply padding to bounds
    minX = minX - settings.CamPadding(1); % L
    maxX = maxX + settings.CamPadding(2); % R
    minY = minY - settings.CamPadding(3); % B
    maxY = maxY + settings.CamPadding(4); % T

    % make the bounds square
    if maxX - minX < maxY - minY
        centerX = (maxX + minX) / 2;
        extent = (maxY - minY) / 2;
        maxX = centerX + extent;
        minX = centerX - extent;
    elseif maxX - minX > maxY - minY
        centerY = (maxY + minY) / 2;
        extent = (maxX - minX) / 2;
        maxY = centerY + extent;
        minY = centerY - extent;
    end

    % prepare the figure
    figure(1);
    mainFig = gcf;
    clf;
    new = ~strcmp(get(mainFig, 'Name'), 'Simulation');
    
    
    % place the window at a convenient location
    if new && numel(settings.InitialWindowPosition) == 4 && ~strcmp(mainFig.WindowStyle, 'docked')
        set(mainFig, 'Position', settings.InitialWindowPosition);
    end
    
    if settings.MaximizeWindow
        set(mainFig, 'Position', get(0, 'Screensize'));
    end

    initialCamera = [minX, maxX, minY, maxY];
    
    ax1 = axes( 'Position', [ 0 0 1 1 ] );
    xlabel('X');
    ylabel('Y');
    axis([minX , maxX , minY , maxY ]);
    axis off;
    axesList = [ax1];

    if settings.PlotEDotHist
        ax2 = axes( 'Position', [ 0.05 0.70 0.3 0.25 ] );
        axes(ax1);
        ax3 = axes( 'Position', [ 0.60 0.70 0.3 0.25 ] );
        axes(ax1);
        axesList = [axesList, ax2,ax3];
    elseif settings.PlotPhiHist || settings.PlotRateWork || settings.PlotApproxRateWork
        ax2 = axes( 'Position', [ 0.05 0.70 0.3 0.25 ] );
        axes(ax1);
        axesList = [axesList, ax2];
    elseif settings.PlotPolarDecomposition
        ax2 = axes( 'Position', [ 0.05 0.70 0.3 0.25 ] );
        axes(ax1); 
        ax3 = axes( 'Position', [ 0.60 0.70 0.3 0.25 ] );
        axes(ax1);
        axesList = [axesList, ax2,ax3];
    elseif settings.PlotContactForceProfile
        axCFP = axes( 'Position', [ 0.15 0.15 0.3 0.25 ] );
        axes(ax1);
        axesList = [axesList, axCFP];
    else 
        a = gca;
        a.Position = [ 0 0 1 1 ];   %%  This is what gives us tight axes in the window!!
    end
    set(mainFig, 'Name', 'Simulation');
    set(gca,'DataAspectRatio',[1 1 1],...
        'PlotBoxAspectRatio',[1 1 1]);

end

