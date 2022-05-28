cla;
clear;
close all;
h = 0.001; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
clear meshes;

%cube
rho = 40;
nu = 0.35;      % Poisson ratio: close to 0.5 for rubber
k = 1e5;     % Young's modulus: 0.01e9 approximate for rubber
[mu,lambda] = toLame(nu,k); 
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.0;  % Rayleigh factor on K
tMaterial1 = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

%matress
rho = 10;
nu = 0.4;      % Poisson ratio: close to 0.5 for rubber
k = 1e3;     % Young's modulus: 0.01e9 approximate for rubber
[mu,lambda] = toLame(nu,k); 
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.1;  % Rayleigh factor on K
tMaterial2 = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

cube = meshLoader("bowlingBall", [], tMaterial1, [0.5,0.5,0.5], false, settings);

s = 4.5;
dist = 8;

baseBeam = meshLoader("mattressHD", [], tMaterial2, [0.07,0.07,0.1], false, settings);
baseBeam.setRigidTransform([90,90,0],[-.1,0,-0.5]);

%bin the bottom
zPos = baseBeam.p(3:3:end);
pinInd = find(zPos <= min(zPos) + 0.1);
baseBeam.pin(sort(pinInd));
%left
yPos = baseBeam.p(2:3:end);
pinInd = find(yPos <= min(yPos) + 0.1);
baseBeam.pin(sort(pinInd));
pinInd = find(yPos >= max(yPos) - 0.1);
baseBeam.pin(sort(pinInd));
%right
xPos = baseBeam.p(1:3:end);
pinInd = find(xPos <= min(xPos) + 0.1);
baseBeam.pin(sort(pinInd));
pinInd = find(xPos >= max(xPos) - 0.1);
baseBeam.pin(sort(pinInd));

cube.setRigidTransform([0,0,0],[0,s,5]);
ncube = Mesh3D(cube);
baseBeam.mergeMesh(ncube);
 
cube.setRigidTransform([0,0,0],[0,-dist,10]);
ncube2 = Mesh3D(cube);
baseBeam.mergeMesh(ncube2);

meshes = AdaptiveMesh3D(baseBeam);


rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 1e-4;
rigidificator.ElastificationThreshold = 1e-3; 
integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.1, 1e-3);
integrator.Gravity = -5;

energyModel = NeoHookean3DEnergy();

meshMeshContactFinder = MeshSCD3D(0.8,true, []);
contactFinder = {meshMeshContactFinder};

settings.MakeVideo = 1;
settings.FramesToRecord = 4000;
settings.PlotSkip = 10;
settings.campos=[12,16,8];
settings.projection = "orthographic";
% settings.PlotEDotHist = 1;
settings.PGSiterations = 50;
settings.SceneName = 'mattressWaking';
settings.OBJDir = './objs/mattressWaking/';
settings.WriteOBJs = true;

td = simulate3D(meshes,h,contactFinder, integrator, rigidificator, settings, energyModel);
writeTDcsv(td,"mattressWaking",["_adaptive"]);