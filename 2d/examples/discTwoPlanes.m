clear 
%close all;
h = 0.01; % time step

rho = 10;
nu = 0.30;
E = 5e4;    
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.00;%3;
alpha1 = 0.002;%1;
material = TriangleMaterial( rho, mu, lambda, alpha0, alpha1 );


%m1 = loadMeshFromPOLY( 'data/fatWiggleRodLongHD.1', material );
%m1.setRigidTransform( 180, [0,0] );

settings = SimulationSettings();
settings.PGSiterations = 500;
settings.WarmStartEnabled = true;
settings.DrawLambdas = true;

settings.CamPadding = [ 1 1 1 1 ]; % L R B T
settings.RunRightAway = false;
settings.SceneName = "discTwoPlanes";
settings.recomputeCacheAinv = true;
settings.warmStartFromQuickSolve= false;
settings.DrawContact = true;
settings.DrawContactFrames = true;

resetMesh = false;
scale = [1,1];
rot = 0;
mesh2d = fetchPoly2D('disc.1',resetMesh, material, scale, rot, settings);
mesh2d.setRigidTransform( 5, [0,0] );
mesh2d.setRigidMotion( 0, [0, 0] ); 
    
mesh2da = AdaptiveMesh( mesh2d );

rigid = EDotMexRigidificator();
rigid.RigidificationThreshold = 1e-4;
rigid.ElastificationThreshold = 5e-3; 
rigid.FrameCount = 3; % number of frames before rigidifying.

integrator = LDLBackwardEuler();
integrator.setComplianceAndBaumgarteFromERPandCFM( h, 0.1, 0.01 );
% integrator.fullElasticContacts = false;
integrator.Gravity = -9.8;

pcf1 = PlaneContactFinder( [0.0, -1], [0,  0.48], 5 ); % normal, position, mu
pcf2 = PlaneContactFinder( [0.0,  1], [0, -0.48], 5 ); % normal, position, mu

% ccf1 = MovingCircleContactFinder( innerRadius, [0, 0], 1 );     % just a bit bigger than the hole
%     ccf1.cfun = @(t) [0,0]; %[ 0, 0.5*(cos(t*8)-1) * (mod( t*8, pi*2*10 ) < pi*2) ];
%     ccf1.dcdt = @(t) [0,0]; %[ 0, -0.5*sin(t*8)*8  * (mod( t*8, pi*2*10 ) < pi*2) ];
%     ccf1.thetafun = @(t)  .5*(sin(t*3 + pi/2) -1)  .* (mod( t*3, pi*2*2 ) < pi*2);
%     ccf1.dthetadt = @(t)  .5*cos(t*3 + pi/2)*3     .* (mod( t*3, pi*2*2 ) < pi*2);

%ccf1 = CircleContactFinder( .15, [-0.4, -0.5], 1 ); % radius, center, mu 
%ccf2 = CircleContactFinder( .15, [ 0.4, -0.5], 1 );  

cfs = { pcf1, pcf2 }%; pcf };
%cfs = { ccf1, ccf2 };

mcf = MeshSCD( 3, true ); % true for quality collisions

% meshes = {mesh2d};
meshes = {mesh2d, mesh2da};  % compare!
td = simulate( meshes, integrator, h, settings, rigid, cfs );

%return

if ( numel(td) > 1 ) 
    figure(4); clf;
    td{1}.plotCompareSiggraphVersion( td{2} ,h);
    pbaspect([3 2 1])
    title('Wall-Clock Simulation Time', 'FontName', 'Linux Biolinum O');
    xlabel('Time(s)', 'FontName', 'Linux Biolinum O');
    ylabel('Time(ms)', 'FontName', 'Linux Biolinum O');
    saveas(gcf,'out/simulation.pdf');
    
%     figure(6); clf;
%     plot(mesha.N-td{2}.logCounts(2,4:end),'Color', [204/255, 85/255, 0], 'LineWidth', 2);
%     title('Number of Rigid Particles per Frame', 'FontName', 'Linux Biolinum O');
%     xlabel('Time(s)', 'FontName', 'Linux Biolinum O');
%     ylabel('number of particles', 'FontName', 'Linux Biolinum O');
%     pbaspect([3 2 1]);
%     saveas(gcf,'out/rigidParticles.pdf');

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