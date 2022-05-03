h = 0.01; % time step

rho = 100;
nu = 0.4; % Poisson ratio: close to 0.5 for rubber
E = 1e5; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.001; % Rayleigh factor on M
alpha1 = 0.01; % Rayleigh factor on K
material = TriangleMaterial( rho, mu, lambda, alpha0, alpha1 );

mesh2d = loadMeshFromPOLY('2d/data/cantileverP05', material)
mesha = AdaptiveMesh(mesh2d);

settings = SimulationSettings();
settings.DrawTimings = 1;
settings.DrawContact = 0;
% settings.DrawLambdas = 1;
settings.CamPadding(3) = 3; % L R B T
settings.DrawEdges = true;
settings.RunRightAway = false;
% settings.PlotEDotHist = 1;
% settings.MakeVideo = true;
%settings.FramesToRecord = 550;
% settings.quicksolveSimulation = 1;
%settings.PCGiterations = 100;
settings.SceneName = "littlePlanetHD";
%settings.PlotPhiHist = 1;
settings.recomputeCacheAinv = true;

rigid = EDotMexRigidificator();
rigid.FrameCount = 3;
rigid.RigidificationThreshold = 5e-4;
rigid.ElastificationThreshold = 5e-3;
rigid.RigidificationThreshold = 1e-5;
rigid.ElastificationThreshold = 1e-4;
rigid.Preconditionner = 'ICHOLICT';
rigid.Permutation = 'DISSECT';

% rigid = EDotMemoryRigidificator();
% rigid.RigidificationThreshold = 1e-5;
% rigid.ElastificationThreshold = 1e-5;
% rigid.ElastificationMult = 1.1;
% rigid.FrameCount = 5; %number of frames before rigidifying.
% rigid.Preconditionner = 'ICHOLICT';
% rigid.Permutation = 'DISSECT';

integrator = LDLBackwardEuler();
ERP = 0.5;
CFM = 0.001;
integrator.setComplianceAndBaumgarteFromERPandCFM( h, ERP, CFM );
integrator.Gravity = -9.8; 

%ccf = CircleContactFinder( 3.0 , [0, -6], 0.3 );
ccf = MovingCircleContactFinder( 3.0 , [0, -6], 0.3 );

    ccf.cfun = @(t) [ 0, 0.3*(cos(t*8)-1) * (mod( t*8, pi*2*10 ) < pi*2) ];
    ccf.dcdt = @(t) [ 0, -0.3*sin(t*8)*8  * (mod( t*8, pi*2*10 ) < pi*2) ];

    ccf.thetafun = @(t)  0;
    ccf.dthetadt = @(t)  0;

td = simulate( { mesha}, integrator, h, settings, {rigid}, ccf );

if ( numel(td) > 1 ) 
    figure(4); clf;
    td{1}.plotCompare( td{2}, 'Rng', 2 );  
    title('Simulation Time');
    figure(5); clf;
    td{1}.plotCompare( td{2});
    title('Simulation Breakdown');
else 
    figure(4); clf;
    td{1}.plot( 'Rng', 2  );  
    title('Simulation Time');
    figure(5); clf;
    td{1}.plot();
    title('Simulation Breakdown');
end