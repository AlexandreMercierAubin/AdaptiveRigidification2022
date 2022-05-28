clear;

h = 0.01; % time step

% MESHES

rho = 10;
nu = 0.38;      % Poisson ratio: close to 0.5 for rubber
E = 7e4;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, E );
alpha0 = 0.001;   % Rayleigh factor on M
alpha1 = 0.3;  % Rayleigh factor on K
tMaterial = TriangleMaterial(rho, mu, lambda, alpha0, alpha1);

mesh = meshLoader("beam", [], tMaterial, [1,1,1], false, settings);
mesh.setRigidTransform( [0,0,0], [2,0,0]); % put it higher in the world frame
% pinning tris
xPos = mesh.p(1:3:end);
pinInd = find(xPos(:,1) == max(xPos));
mesh.pin(sort(pinInd));

%mesh.setRigidMotion( [0,0.05,1], [0,0,0] );

mesha = AdaptiveMesh3D(mesh);

% rigidificator = EDot3DMexRigidificator();
rigidificator = EDot3DRigidificator();
rigidificator.RigidificationThreshold = 1e-8;
rigidificator.ElastificationThreshold = 1e-7; 
rigidificator.FrameCount = 3;

integrator = LDLBackwardEuler3D();
integrator.Compliance = 0.0;

energyModel = NeoHookean3DEnergy();

planeContactFinder = NullContactFinder(3);
contactFinder = {planeContactFinder};

% CONFIG

settings = Simulation3DSettings();

% settings.PlotEDotHist = 1;

settings.PlotSkip = 5;
settings.campos=[-30,40,5]*.5;

settings.camtarget=[1,1,0];

% settings.MakeVideo = 1;
settings.SceneName = 'pinnedbeam3d-a';
settings.FramesToRecord = 2000;
settings.WriteOBJs = 1;

td = simulate3D( mesha, h, contactFinder, integrator, rigidificator, settings, energyModel);
