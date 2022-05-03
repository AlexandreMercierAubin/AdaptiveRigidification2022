clear;

h = 0.1; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

rho = 10;
nu = 0.33;% 0.33;      % Poisson ratio: close to 0.5 for rubber
E = 1e5;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, E );
alpha0 = 0;%.0001;   % Rayleigh factor on M
alpha1 = 0.01;  % Rayleigh factor on K
tMaterial = TriangleMaterial(rho, mu, lambda, alpha0, alpha1);

mesh = meshLoader("beam", [], tMaterial, [1,1,1], false, settings);
mesh.setRigidTransform( [0,0,0], [0,0,1.25]); % put it higher in the world frame
% pinning tris
xPos = mesh.p(1:3:end);
pinInd = find(xPos(:,1) == max(xPos));
%mesh.pin(sort(pinInd));

mesh.setRigidMotion( [0,0.5,1], [0,0,0] );

mesha = AdaptiveMesh3D(mesh);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 1e-4;
rigidificator.ElastificationThreshold = 5e-3;% 5e-1; 
rigidificator.FrameCount = 5;

integrator = BackwardEuler3D();
integrator.Compliance = 0.0;
integrator.Gravity = 0;

energyModel = StVenantKirchoff3DEnergy();

planeContactFinder = NullContactFinder(3);
contactFinder = {planeContactFinder};

settings.WriteOBJs = 0;
% settings.DrawRigidDOFs = 1;
% settings.PlotSkip = 2;
settings.campos=[-17,20,5];
% settings.PlotEDotHist = 1;
settings.MakeVideo = 1;
settings.SceneName = 'spinBeam';
settings.FramesToRecord = 500;
% settings.DrawForces = true;
simulate3D( mesha, h, contactFinder, integrator, rigidificator, settings, energyModel);
