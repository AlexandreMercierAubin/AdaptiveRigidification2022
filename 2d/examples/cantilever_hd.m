clear

gravity = -9.8; % gravity
h = 0.01; % time step

rho = 10;
nu = 0.33; % Poisson ratio: close to 0.5 for rubber
E = 6e3; % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, E );

% Set damping manually
alpha0 = 0.001; % Rayleigh factor on M
alpha1 = 0.0; % Rayleigh factor on K

% Set damping from damping factor
% compute first modal frequency via eigs
% where omega0 is sqrt(D(1,1))
% (set a breakpoint in LDLBackwardEuler to harvest this)
% Mii=M(ii,ii);Kii=(K(ii,ii)+K(ii,ii)')/2;
% [V,D] = eigs(Kii,Mii,7,'smallestabs'); 
% omeag0 = sqrt(-D(1,1))
% If unpinned 2D then use D(4,4)
% If unpinned 3D then use D(7,7)
omega0 = 6.7485; % divide by 2pi to get period
df = 0.2;  % zero is undamped, 1 is critical
% set alpha1 from alpha0
alpha0 = 0;
alpha1 = (2*df-alpha0/omega0)/omega0;
% set alpha0 from alpha1
% alpha1 = 0;
% alpha0 = (2*df-alpha1*omega0)*omega0;


material = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

settings = SimulationSettings();

% settings.MakeVideo = 1;
settings.SceneName = 'cantilever_hd';
% settings.DrawTimings = 1;
% settings.PrintTimings = 0;
settings.CamPadding(3) = 3;
% settings.plotTriImplicitRigidificationElastification = 1;
settings.PlotEDotHist = 1;
settings.RecordFramePositionInTD = true;
% settings.DrawDv = true;
% settings.DrawApproxDv = true;

resetMesh = false;
scale = [1,1];
rot = 0;
mesh2d = fetchPoly2D('cantileverP2',resetMesh, material, scale, rot, settings);

% pin the left hand side
px = mesh2d.p(1:2:end);
minx = min( px );
mesh2d.pin( find( px < minx + 0.1 ) );

mesha = AdaptiveMesh( mesh2d );

rigid = EDotMexRigidificator();
rigid.RigidificationThreshold = 1e-6;
rigid.ElastificationThreshold = 1e-5;

integrator = LDLBackwardEuler();
integrator.Gravity = gravity*0.5;
integrator.separateQuicksolveGravity = false;

settings.PlotSkip = 3;

td = simulate( {mesh2d, mesha}, integrator, h, settings, rigid );
