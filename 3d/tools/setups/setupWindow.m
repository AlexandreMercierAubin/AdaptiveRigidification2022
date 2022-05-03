function [mainFig, axesList, initialCamera] = setupWindow(settings,meshes)
%[mainFig, axesList, initialCamera] = setupWindow(settings,meshes)
%sets up the figure window to plot the scene with respect to the settings.
minX = 0;
minY = 0;
minZ = 0;
maxX = 0;
maxY = 0;
maxZ = 0;

for k = 1:size(meshes, 2)
    minX = min(minX, min(meshes{k}.p(1:3:end)));
    minY = min(minY, min(meshes{k}.p(2:3:end)));
    minZ = min(minZ, min(meshes{k}.p(3:3:end)));
    maxX = max(maxX, max(meshes{k}.p(1:3:end)));
    maxY = max(maxY, max(meshes{k}.p(2:3:end)));
    maxZ = max(maxZ, max(meshes{k}.p(3:3:end)));
end

extent = max((maxX - minX) / 2, (maxY - minY) / 2);
extent = max(extent, (maxZ - minZ) / 2);

centerX = (maxX + minX) / 2;
centerY = (maxY + minY) / 2;
centerZ = (maxZ + minZ) / 2;

maxX = centerX + extent;
minX = centerX - extent;
maxY = centerY + extent;
minY = centerY - extent;
maxZ = centerZ + extent;
minZ = centerZ - extent;

% prepare the figure
figure(1);
mainFig = gcf;
set(mainFig, 'ToolBar', 'none');
clf;
new = ~strcmp(get(mainFig, 'Name'), 'Simulation');

% place the window at a convenient location
if new && numel(settings.InitialWindowPosition) == 4 && ~strcmp(mainFig.WindowStyle, 'docked')
    set(mainFig, 'Position', settings.InitialWindowPosition);
end    

disp("-----------------------------------------");
disp("Keyboard Controls                        ");
disp("-----------------------------------------");
disp("r               reset                    ");
disp("p space         start and stop simulation");
disp("s               step once                ");
disp("e               toggle elastification    ");
disp("d               toggle rigidification    ");
disp("escape          quit simulation          ");
disp("-----------------------------------------");


%axes;
ax1 = axes('Units', 'normalized', 'Position', [ 0 0 1 1 ] );
camtarget(ax1, settings.camtarget );
camproj(ax1, settings.projection);
ax1.CameraViewAngle = settings.camfov;
ax1.CameraViewAngleMode = 'manual';
disableDefaultInteractivity(ax1);
campos(ax1, settings.campos);
if settings.Shading
    camlight(ax1,'right', 'infinite');
    shading(ax1,'interp'); % interp;
    lighting(ax1,'gouraud'); % gouraud;
end
%     axis vis3d equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
initialCamera = [minX - 5, maxX + 5, minY - 5, maxY + 5, minZ - 5, maxZ + 5];
axis(ax1,initialCamera);
axis off;
ax1.Clipping = 'off'; 
cameratoolbar('Show');
%     hold( ax1, 'on' );
axesList = [ax1];

if ( settings.PlotEDotHist )
    ax2 = axes( 'Position', [ 0.05 0.70 0.3 0.25 ] );
    disableDefaultInteractivity(ax2);
    axes(ax1);
    ax5 = axes( 'Position', [ 0.60 0.70 0.3 0.25 ] );
    disableDefaultInteractivity(ax5);
    axes(ax1);
    axesList=[axesList,ax2,ax5];
elseif ( settings.PlotMomentum )
    ax3 = axes( 'Position', [ 0.05 (.5-0.25*.5) 0.3 0.25 ] );
    disableDefaultInteractivity(ax3);
    ax4 = axes( 'Position', [ 0.05 0.05 0.3 0.25 ] );
    disableDefaultInteractivity(ax4);
    axes(ax1);
    axesList=[axesList,ax3,ax4];
end
set(mainFig, 'Name', 'Simulation');
end

