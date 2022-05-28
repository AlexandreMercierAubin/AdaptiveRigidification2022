cla;
clear;
close all;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
clear meshes;

%cube
rho = 5;
nu = 0.35;      % Poisson ratio: close to 0.5 for rubber
k = 5e3;     % Young's modulus: 0.01e9 approximate for rubber
[mu,lambda] = toLame(nu,k); 
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.0;  % Rayleigh factor on K
tMaterial1 = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

%matress
rho = 5;
nu = 0.4;      % Poisson ratio: close to 0.5 for rubber
k = 1.5e3;     % Young's modulus: 0.01e9 approximate for rubber
[mu,lambda] = toLame(nu,k); 
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.06;  % Rayleigh factor on K
tMaterial2 = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

cube = meshLoader("beam28", [], tMaterial1, [0.15,0.15,0.05], false, settings);

s = 6;
dist = 3;
heightChange = 3;

baseBeam = meshLoader("beam913", [], tMaterial2, [0.6,0.6,0.6]);
baseBeam.setRigidTransform([0,90,0],[-.1,0,0]);

cube.setRigidTransform([0,0,0],[0,0,2]);
baseMesh = Mesh3D(baseBeam);
baseMesh.mergeMesh(cube);
baseMesh.setRigidTransform([0,0,0],[0,s,1]);
zPos = baseMesh.p(3:3:end);
pinInd = find(zPos == min(zPos));
baseMesh.pin(sort(pinInd));
meshes1 = AdaptiveMesh3D(baseMesh);

cube.setRigidTransform([0,0,0],[0,0,heightChange]);
baseMesh = Mesh3D(baseBeam);
baseMesh.mergeMesh(cube);
baseMesh.setRigidTransform([0,0,0],[0,s-dist,1]);
baseMesh.pin(sort(pinInd));
meshes2 = AdaptiveMesh3D(baseMesh);

cube.setRigidTransform([0,0,0],[0,0,heightChange]);
baseMesh = Mesh3D(baseBeam);
baseMesh.mergeMesh(cube);
baseMesh.setRigidTransform([0,0,0],[0,s-dist*2,1]);
baseMesh.pin(sort(pinInd));
meshes3 = AdaptiveMesh3D(baseMesh);


rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 7e-6;
rigidificator.ElastificationThreshold = 7e-5; 
integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.05, 1e-3);
integrator.Gravity = -8;

energyModel = NeoHookean3DEnergy();

meshMeshContactFinder = MeshSCD3D(0.8,true, []);
contactFinder = {meshMeshContactFinder};

settings.MakeVideo = 1;
settings.FramesToRecord = 250;
% settings.PGSiterations = 20;
% settings.DrawVelocities = true;
% settings.DrawForces = true;
% settings.DrawLambdas = 1;
% settings.DrawRigidDOFs = 1;
settings.campos=[12,16,8];
settings.projection = "orthographic";
% settings.PlotEDotHist = 1;
settings.PGSiterations = 50;
settings.SceneName = 'beamWaking';
settings.OBJDir = './objs/beamWaking/';
settings.WriteOBJs = true;
meshes = {meshes1,meshes2,meshes3};
simulate3D(meshes,h,contactFinder, integrator, rigidificator, settings, energyModel);