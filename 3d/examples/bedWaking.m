cla;
clear;
close all;
h = 0.001; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
clear meshes;

%cube
rho = 10;
nu = 0.35;      % Poisson ratio: close to 0.5 for rubber
k = 5e3;     % Young's modulus: 0.01e9 approximate for rubber
[mu,lambda] = toLame(nu,k); 
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.0;  % Rayleigh factor on K
tMaterial1 = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

%matress
rho = 10;
nu = 0.4;      % Poisson ratio: close to 0.5 for rubber
k = 1.5e3;     % Young's modulus: 0.01e9 approximate for rubber
[mu,lambda] = toLame(nu,k); 
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.06;  % Rayleigh factor on K
tMaterial2 = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

cube = meshLoader("bowlingBall", [], tMaterial1, [0.15,0.15,0.15], false, settings);

s = 6;
dist = 3;
heightChange = 4;

baseBeam = meshLoader("bed", [], tMaterial2, [0.8,0.8,1]);
baseBeam.setRigidTransform([90,90,0],[-.1,0,0]);

cube.setRigidTransform([0,0,0],[0,0,2.5]);
baseMesh = Mesh3D(baseBeam);
baseMesh.mergeMesh(cube);
baseMesh.setRigidTransform([0,0,0],[0,s,1]);

%bin the bottom
zPos = baseMesh.p(3:3:end);
pinInd = find(zPos <= min(zPos) + 0.05);
baseMesh.pin(sort(pinInd));
%pin the wooden part
yPos = baseMesh.p(1:3:end);
pinInd2 = find(yPos <= min(yPos) + 0.05);
baseMesh.pin(sort(pinInd2));

meshes1 = AdaptiveMesh3D(baseMesh);

cube.setRigidTransform([0,0,0],[0,0,heightChange]);
baseMesh = Mesh3D(baseBeam);
baseMesh.mergeMesh(cube);
baseMesh.setRigidTransform([0,0,0],[0,s-dist,1]);
baseMesh.pin(sort(pinInd));
baseMesh.pin(sort(pinInd2));
meshes2 = AdaptiveMesh3D(baseMesh);

cube.setRigidTransform([0,0,0],[0,0,heightChange]);
baseMesh = Mesh3D(baseBeam);
baseMesh.mergeMesh(cube);
baseMesh.setRigidTransform([0,0,0],[0,s-dist*2,1]);
baseMesh.pin(sort(pinInd));
baseMesh.pin(sort(pinInd2));
meshes3 = AdaptiveMesh3D(baseMesh);


rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 7e-6;
rigidificator.ElastificationThreshold = 1e-4; 
integrator = BackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.1, 1e-3);
integrator.Gravity = -8;

energyModel = NeoHookean3DEnergy();

meshMeshContactFinder = MeshSCD3D(0.8,true, []);
contactFinder = {meshMeshContactFinder};

settings.MakeVideo = 1;
settings.FramesToRecord = 500;
% settings.PGSiterations = 20;
% settings.DrawVelocities = true;
% settings.DrawForces = true;
% settings.DrawLambdas = 1;
% settings.DrawRigidDOFs = 1;
settings.campos=[12,16,8];
settings.projection = "orthographic";
% settings.PlotEDotHist = 1;
settings.PGSiterations = 50;
settings.SceneName = 'bedWaking';
settings.OBJDir = './objs/bedWaking/';
settings.WriteOBJs = true;
meshes = {meshes1,meshes2,meshes3};
simulate3D(meshes,h,contactFinder, integrator, rigidificator, settings, energyModel);