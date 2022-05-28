cla;
clear;
close all;
h = 0.002; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear meshes;

rho = 50;
nu = 0.35;      % Poisson ratio: close to 0.5 for rubber
k = 5e3;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu,k );
alpha0 = 0.001;   % Rayleigh factor on M
alpha1 = 0.001;  % Rayleigh factor on K
tMaterialArmadillo = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];
tMaterialArmadillo.cacheName = 'armadilloSingleMat';

% settings.MakeVideo = 1;
settings.SceneName = 'singleArmadillo';
% settings.FramesToRecord = 500;
% settings.InitialWindowPosition = [0,0,1920,1080];
settings.campos=[5,10,3];

resetMesh = false;
baseMesh = polyLoader("arma_6", tMaterialArmadillo,[1,1,1],resetMesh, settings);

baseMesh.setRigidTransform([90,0,0],[0,0,0.4]);

meshes = AdaptiveMesh3D(baseMesh);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 5e-5;
rigidificator.ElastificationThreshold = 5e-4; 
integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0.1, 0.1 );
integrator.Gravity = -9.8;

planeContactFinder = PlaneContactFinder3D([0,0,1], [0,0,0], 0.7);
meshMeshContactFinder = MeshSCD3D(0.0,false, []);
contactFinder = {planeContactFinder, meshMeshContactFinder};

energyModel = NeoHookean3DEnergy();

simulate3D(meshes,h,contactFinder, integrator, rigidificator, settings, energyModel);