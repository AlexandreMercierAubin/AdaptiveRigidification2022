clear all;
% close all;

h = 0.01; % time step

rho = 1;
nu = 0.4; % Poisson ratio: close to 0.5 for rubber
E = 5e3; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.001; % Rayleigh factor on M
alpha1 = 0.01; % Rayleigh factor on K
material = TriangleMaterial( rho, mu, lambda, alpha0, alpha1 );

settings = SimulationSettings();
settings.DrawTimings = 1;
settings.CamPadding(3) = 3;
settings.SceneName = 'ground_cantilever_mega';
%settings.MakeVideo = true;
%settings.DrawEdges = true;
%settings.DrawLambdas = 1;

% settings.quicksolveSimulation = true;
settings.PCGiterations = 1;
settings.PGSiterations = 100;
settings.recomputeCacheAinv = true;

resetMesh = false;
scale = [1,1];
rot = 0;
mesh2d = AdaptiveMesh(fetchPoly2D('cantileverP05',resetMesh, material, scale, rot, settings));

rigid = EDotMexRigidificator();
rigid.RigidificationThreshold = 1e-6;
rigid.ElastificationThreshold = 1e-5;

integrator = LDLBackwardEuler();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.1, 0.001 );
integrator.Gravity = -9.8; 
% integrator.useFullAinv = true;

pcf = PlaneContactFinder( [0.0, 1.0] , [0, -3], 0.0 );

td = simulate( mesh2d, integrator, h, settings, rigid, pcf );


if ( numel(mesh2d) >1 ) 
    figure(4); clf;
    td{1}.plotCompare( td{2}, 'Rng', 2 );  
    title('Simulation Time');
    figure(5); clf;
    td{1}.plotCompare( td{2} );
    title('Simulation Breakdown');
else 
    figure(4); clf;
    td{1}.plot( 'Rng', 2  );  
    title('Simulation Time');
    figure(5); clf;
    td{1}.plot();
    title('Simulation Breakdown');
end

