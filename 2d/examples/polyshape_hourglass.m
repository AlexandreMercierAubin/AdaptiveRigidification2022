clear;
close all;
h = 0.01; % time step

rho = 100;
nu = 0.35; % Poisson ratio: close to 0.5 for rubber
E = 5e4; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.001; % Rayleigh factor on M
alpha1 = 0.1; % Rayleigh factor on K
material = TriangleMaterial( rho, mu, lambda, alpha0, alpha1 );

settings = SimulationSettings();
settings.CamPadding(3) = 3; % L R B T
settings.MakeVideo = 1;
settings.SceneName = 'hourglass';
settings.FramesToRecord = 500;
settings.PGSiterations = 100;

resetMesh = false;
scale = [1,1];
rot = 90;
mesh2d = AdaptiveMesh(fetchPoly2D('cantileverP3',resetMesh, material, scale, rot, settings));
mesh2d.setRigidTransform(0,[-1.5,2]);
mesha = AdaptiveMesh( mesh2d );

% settings.PlotEDotHist = true;
% settings.InitialWindowPosition = [0,0,1920,1080];

rigid = EDotMexRigidificator();
rigid.RigidificationThreshold = 1e-3;
rigid.ElastificationThreshold = 1e-2;

integrator = LDLBackwardEuler();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.1,0.001);
integrator.Gravity = -9.8; 

friction = 0;

%nodes must only be on edges and be ordered.
[I,J] = readNODE('2d/data/cup.node');
I(:,1) = I(:,1)*2;
rotation = 90;
R = [cosd(rotation), - sind(rotation);...
         sind(rotation), cosd(rotation)];
leftI = I*R';
leftI(:,1) = leftI(:,1) +0.1;
leftI(:,2) = leftI(:,2) -4;

rotation = 270;
R = [cosd(rotation), - sind(rotation);...
         sind(rotation), cosd(rotation)];
rightI = I*R';
rightI(:,1) = rightI(:,1) + 0.9;
rightI(:,2) = rightI(:,2) - 4;

contactFinder = PolyShapeContactFinder( leftI, friction );
contactFinder2 = PolyShapeContactFinder( rightI, friction );

settings.MakeVideo = 1;
settings.CamPadding = [1,1,10,1];

td = simulate( {mesh2d,mesha}, integrator, h, settings, rigid, {contactFinder,contactFinder2} );
% 
% if ( numel(td) > 1 ) 
%     figure(4); clf;
%     td{1}.plotCompareSiggraphVersion( td{2} ,h);
%     pbaspect([3 2 1])
%     title('Wall-Clock Simulation Time', 'FontName', 'Linux Biolinum O');
%     xlabel('Time(s)', 'FontName', 'Linux Biolinum O');
%     ylabel('Time(ms)', 'FontName', 'Linux Biolinum O');
%     saveas(gcf,'out/simulation.pdf');
% 
%     figure(7); clf;
%     hold on;
%     plot(td{2}.logCounts(7,4:end),'Color', [204/255, 85/255, 0], 'LineWidth', 2 );
% %     plot(td{1}.logCounts(7,4:end),'Color', [170/255, 20/255, 20/255], 'LineWidth', 2 );
%     title('Variation in Degrees of Freedom Over Time', 'FontName', 'Linux Biolinum O');
%     xlabel('Time(s)', 'FontName', 'Linux Biolinum O');
%     ylabel('DOFs', 'FontName', 'Linux Biolinum O');
%     pbaspect([3 2 1]);
%     saveas(gcf,'out/totalDofs.pdf');
% end
