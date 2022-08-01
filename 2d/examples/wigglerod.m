clear 
close all;
h = 0.01; % time step

settings = SimulationSettings();
settings.PGSiterations = 20;
settings.CamPadding = [ 0 0 0 0 ]; % L R B T
% settings.DrawTimings = true;
settings.RunRightAway = true;
%settings.MakeVideo = 1;
%settings.FramesToRecord = 500;
settings.SceneName = "wigglerodLongHD";
settings.WarmStartEnabled = true;

rho = 30;
nu = 0.30;
E = 1e5;    
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.00;%3;
alpha1 = 0.1;
material = TriangleMaterial( rho, mu, lambda, alpha0, alpha1 );


%m1 = loadMeshFromPOLY( 'data/fatWiggleRodLongHD.1', material );
%m1.setRigidTransform( 180, [0,0] );
resetMesh = false;
scale = [1,1];
rot = 0;
m1 = fetchPoly2D('wigglerodlongHD.1',resetMesh, material, scale, rot, settings);
 % this was probably not supposed to be adaptive for the comparison!
m1a = AdaptiveMesh( m1 );
 
rigid = EDotMexRigidificator();
rigid.RigidificationThreshold = 1e-6;
rigid.ElastificationThreshold = 2e-5; 
rigid.FrameCount = 3; % number of frames before rigidifying.

integrator = LDLBackwardEuler();
ERP = 0.2;
CFM = 0.001;
integrator.setComplianceAndBaumgarteFromERPandCFM( h, ERP , CFM );

integrator.Gravity = -9.8;

%pcf = PlaneContactFinder( [0.0, 1] , [0, -0.5], 1 ); % normal, position, mu
pcf = PlaneContactFinder( [0.0, 1] , [0, -0.8], 1 ); % normal, position, mu

% ccf1 = MovingCircleContactFinder( innerRadius, [0, 0], 1 );     % just a bit bigger than the hole
%     ccf1.cfun = @(t) [0,0]; %[ 0, 0.5*(cos(t*8)-1) * (mod( t*8, pi*2*10 ) < pi*2) ];
%     ccf1.dcdt = @(t) [0,0]; %[ 0, -0.5*sin(t*8)*8  * (mod( t*8, pi*2*10 ) < pi*2) ];
%     ccf1.thetafun = @(t)  .5*(sin(t*3 + pi/2) -1)  .* (mod( t*3, pi*2*2 ) < pi*2);
%     ccf1.dthetadt = @(t)  .5*cos(t*3 + pi/2)*3     .* (mod( t*3, pi*2*2 ) < pi*2);

ccf1 = CircleContactFinder( .15, [-0.4, -0.5], 1 ); % radius, center, mu 
ccf2 = CircleContactFinder( .15, [ 0.4, -0.5], 1 );  

cfs = { pcf, ccf1, ccf2 };
%cfs = { ccf1, ccf2 };

meshes = {m1, m1a};
%meshes = {m1, m1a};  % compare!
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