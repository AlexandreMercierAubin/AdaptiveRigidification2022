clear
close all;
g = -9.8; % gravity
h = 0.01; % time step
friction = 1.0;

rho = 30;
nu = 0.4; % Poisson ratio: close to 0.5 for rubber
E = 2e5; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.001; % Rayleigh factor on M
alpha1 = 0.01;  % Rayleigh factor on K
material1 = TriangleMaterial( rho, mu, lambda, alpha0, alpha1 );

rho = 10;
nu = 0.30; % Poisson ratio: close to 0.5 for rubber
E = 1e5; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.00;%3;
alpha1 = 0.00;%2;
material2 = TriangleMaterial( rho, mu, lambda, alpha0, alpha1, [0.5,0.65,0.5] );

rho = 20;
nu = 0.30; % Poisson ratio: close to 0.5 for rubber
E = 5e3; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.00;%3;
alpha1 = 0.02;%2;
material3 = TriangleMaterial( rho, mu, lambda, alpha0, alpha1, [0.5,0.5,0.65] );

settings = SimulationSettings();
% settings.MakeVideo = 1;
% settings.FramesToRecord = 600;
settings.DrawTimings = 0;
settings.CamPadding(1:2) = [1,1]; % L R B T
settings.DrawContact = 0;
settings.DrawLambdas = 0;
settings.RunRightAway = 1;
settings.PGSiterations = 100;
% settings.InitialWindowPosition = [0,0,1920,1080];
settings.ElasticContacts = 1;  % need to remove this variable I think
% settings.PlotSpyA = 1;
%settings.MakeVideo = 1;
settings.SceneName = 'multimergeTest';
%settings.FramesToRecord = 600;

mesh2d = fetchMesh2D('multimerge_test',false, @multimeshmerge_test_builder, [material1,material2,material3],settings);

% MUST NOT merge adaptive meshes... or write code to do the merge
% properly... first merge then make adaptive.
mesha = AdaptiveMesh( mesh2d );

rigid = EDotMexRigidificator();
rigid.FrameCount = 3;
rigid.RigidificationThreshold = 5e-4;
rigid.ElastificationThreshold = 1e-2;

integrator = LDLBackwardEuler();
ERP = 0.5;
CFM = 0.01;
integrator.setComplianceAndBaumgarteFromERPandCFM( h, ERP, CFM );
integrator.Gravity = g;

contactFinder = PlaneContactFinder([0.0, 1.0] , [0, -1], friction );

mcf = MeshSCD( 1, true ); % true for quality collisions

%meshes = {mesh, mesha};
meshes = mesha;
td = simulate( meshes, integrator, h, settings, rigid, contactFinder, mcf );

return;

if ( numel(td) > 1 ) 
    figure(4); clf;
    td{1}.plotCompareSiggraphVersion( td{2} ,h);
    pbaspect([3 2 1])
    title('Wall-Clock Simulation Time', 'FontName', 'Linux Biolinum O');
    xlabel('Time(s)', 'FontName', 'Linux Biolinum O');
    ylabel('Time(ms)', 'FontName', 'Linux Biolinum O');
    saveas(gcf,'out/simulation.pdf');
   

    figure(7); clf;
    hold on;
    plot(td{2}.logCounts(7,4:end),'Color', [204/255, 85/255, 0], 'LineWidth', 2 );
%     plot(td{1}.logCounts(7,4:end),'Color', [170/255, 20/255, 20/255], 'LineWidth', 2 );
    title('Variation in Degrees of Freedom Over Time', 'FontName', 'Linux Biolinum O');
    xlabel('Time(s)', 'FontName', 'Linux Biolinum O');
    ylabel('DOFs', 'FontName', 'Linux Biolinum O');
    pbaspect([3 2 1]);
    saveas(gcf,'out/totalDofs.pdf');
end