clear 
% close all;
h = 0.01; % time step

rho = 100;
nu = 0.30; % Poisson ratio: close to 0.5 for rubber
k = 8e4; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, k );

alpha0 = 0.0003; % a larger alpha 0 will fight against contact stabilizaiton :(
alpha1 = 0.03;
%alpha0 = 0;
%alpha1 = 0;

%strainUpperBound = 1.2;
%strainLowerBound = 0.8;
material = TriangleMaterial( rho, mu, lambda, alpha0, alpha1);

settings = SimulationSettings();
% settings.MakeVideo = 1;
%settings.FramesToRecord = 200;
settings.DrawTimings = 0;
settings.DrawLambdas = 0;
settings.PGSiterations = 20;
settings.RunRightAway = true;
settings.CamPadding = [0.5 0.5 0 -1.75 ];
settings.PlotEDotHist = true;
% settings.InitialWindowPosition = [0,0,1920,1080];
settings.SceneName = 'IShape';
%settings.DrawApproxEDots = true;
%settings.DrawApproxDv = true;
%settings.PlotPhiHist = 1;
%settings.PlotRateWork = 1;
% settings.PlotApproxRateWork = 1;
settings.recomputeCacheAinv = true;

resetMesh = false;
scale = [1,1];
rot = 0;
%m1 = fetchPoly2D('IShape.1',false, material, scale, rot, settings);
m1 = fetchPoly2D('IShapeHD.1',true, material, scale, rot, settings);
%m1 = AdaptiveMesh( m1 );
m1.setRigidTransform( 0, [-1.5 0] );
m2 = AdaptiveMesh( m1 );
m2.setRigidTransform( 0, [1.5 0] );

% rigid = EDotMexRigidificator();
% rigid.RigidificationThreshold = 1e-5;
% rigid.ElastificationThreshold = 1e-4;
% rigid.Permutation = 'DISSECT';
% rigid.FrameCount = 5; %number of frames before rigidifying.
% rigid.Preconditionner = 'ICHOLICT';

rigid = EDotMexRigidificator();
rigid.RigidificationThreshold = 1e-5;
rigid.ElastificationThreshold = 1e-4;
rigid.Permutation = 'DISSECT';
rigid.FrameCount = 3; %number of frames before rigidifying.
rigid.Preconditionner = 'ICHOLICT';


% Note clear that the approx rigificator (or the memory rigidificator) can
% be made to work in short order.  Keep with the vanilla option above.
integrator = LDLBackwardEuler();
ERP = 0.5;
CFM = 0.1;
integrator.setComplianceAndBaumgarteFromERPandCFM( h, ERP, CFM );
integrator.Gravity = -9.8; 
integrator.separateQuicksolveGravity = true;
% integrator.useFullAinv = true;

frictionCoefficient = 2.0; 
pcf = PlaneContactFinder([0.0, 1.0] , [0.01, -1], frictionCoefficient );

%meshes = {m1,m2};    % run as comparison
meshes = m2;
td = simulate( meshes, integrator, h, settings, rigid, pcf );

return

if ( numel(meshes) >1 ) 
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
