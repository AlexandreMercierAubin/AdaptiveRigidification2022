clear 
%close all;
h = 0.01; % time step

rho = 10;
nu = 0.30;
E = 2e6;    
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.00;%3;
alpha1 = 0.002;%1;
material1 = TriangleMaterial( rho, mu, lambda, alpha0, alpha1 );

rho = 10;
nu = 0.30;
E = 2e4;    
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.00;%3;
alpha1 = 0.002;%1;
material2 = TriangleMaterial( rho, mu, lambda, alpha0, alpha1 );

settings = SimulationSettings();
settings.PGSiterations = 500;
settings.WarmStartEnabled = true;
settings.DrawLambdas = true;
settings.ElasticContacts = 1;   % why is this not on by default!?
settings.CamPadding = [ 1 1 1 1 ]; % L R B T
% settings.DrawTimings = true;
settings.RunRightAway = true;
settings.MakeVideo = 1;
%settings.FramesToRecord = 500;
settings.SceneName = "foldingCVertical.mp4";
%settings.PlotEDotHist = true;
settings.recomputeCacheAinv = true;
settings.warmStartFromQuickSolve= false;
% settings.DrawContactFrames = 1;
settings.DrawContact = true;

mesh2d = fetchMesh2D('foldingCmediumsmall3',false, @foldingVerticalCnodiscBuilder, [material1,material2], settings);

% Editing these scripts is making me think about this:
% https://www.mathworks.com/help/matlab/json-format.html

mesh2da = AdaptiveMesh( mesh2d );

rigid = EDotMexRigidificator();
rigid.RigidificationThreshold = 1e-4;
rigid.ElastificationThreshold = 5e-3; 
rigid.FrameCount = 3; % number of frames before rigidifying.

integrator = LDLBackwardEuler();
integrator.setComplianceAndBaumgarteFromERPandCFM( h, 0.1, 0.01 );
% integrator.fullElasticContacts = false;
integrator.Gravity = -9.8;

cfs = {}%; pcf };
%cfs = { ccf1, ccf2 };

mcf = MeshSCD( 1, true ); % true for quality collisions

meshes = {mesh2da};
td = simulate( meshes, integrator, h, settings, rigid, cfs, mcf );

%return

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