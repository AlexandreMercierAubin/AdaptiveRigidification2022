clear ;
%close all;
h = 0.01; % time step

rho = 100;
nu = 0.30; % Poisson ratio: close to 0.5 for rubber
k = 8e4; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, k );
alpha0 = 0.0003; % a larger alpha 0 will fight against contact stabilizaiton :(
alpha1 = 0.1;
strainUpperBound = 1.2;
strainLowerBound = 0.8;
materialColor = [0.5,0.5,0.5];
material = TriangleMaterial( rho, mu, lambda, alpha0, alpha1,materialColor, strainUpperBound, strainLowerBound);

settings = SimulationSettings();
% settings.MakeVideo = 1;
% settings.FramesToRecord = 200;
settings.DrawTimings = 0;
settings.DrawLambdas = 1;
settings.PGSiterations = 100;
% settings.RunRightAway = false;
settings.CamPadding = [1 1 0 -1.75 ];
% settings.PlotEDotHist = true;
% settings.PlotPhiHist = true;
% settings.InitialWindowPosition = [0,0,1920,1080];
settings.SceneName = 'LShape';
settings.recomputeCacheAinv = true;
% settings.DrawApproxDv = true;


resetMesh = false;
scale = [1,1];
rot = 0;
m1 = fetchPoly2D('LShape.1',true, material, scale, rot, settings);
%m1 = AdaptiveMesh( m1 );
m1.setRigidTransform( 0, [-5 2] );
m2 = AdaptiveMesh( m1 );
m2.setRigidTransform( 0, [5 0] );

rigid = EDotMexRigidificator();
rigid.RigidificationThreshold = 1e-6;
rigid.ElastificationThreshold = 5e-5;

rigid.FrameCount = 3; %number of frames before rigidifying.

integrator = LDLBackwardEuler();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.2, 0.01 );
integrator.Gravity = -9.8; 

pcf = PlaneContactFinder([0.0, 1.0] , [0.01, -1], 1.0);

meshes = {m2};    % run as comparison
td = simulate( meshes, integrator, h, settings, rigid, pcf );