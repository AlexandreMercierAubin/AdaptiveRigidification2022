cla;
clear;
close all;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
clear meshes;
rho = 10;
nu = 0.4;      % Poisson ratio: close to 0.5 for rubber
k = 1e4;     % Young's modulus: 0.01e9 approximate for rubber
[mu,lambda] = toLame(nu,k); 
alpha0 = 0.01;   % Rayleigh factor on M
alpha1 = 0.1;  % Rayleigh factor on K
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

settings.MakeVideo = 1;
settings.FramesToRecord = 200;
settings.InitialWindowPosition = [0,0,1920,1080];
% settings.PGSiterations = 20;
% settings.DrawVelocities = true;
% settings.DrawForces = true;
settings.DrawLambdas = 1;
settings.DrawRigidDOFs = 1;
settings.campos=[10,10,1];
settings.PlotEDotHist = 1;
settings.PGSiterations = 100;
settings.SceneName = 'beamLowDef';

resetMesh = false;

baseMesh = meshLoader("beam28", [], tMaterial, [0.6,0.6,0.6], resetMesh, settings);
baseMesh.setRigidTransform([0,90,0],[0,0,2]);
baseMesh2 = meshLoader("beam28", [], tMaterial, [0.5,0.5,0.5], resetMesh, settings);
baseMesh2.setRigidTransform([0,90,0],[0,0,6]);
baseMesh.mergeMesh(baseMesh2);
meshes = AdaptiveMesh3D(baseMesh);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 1e-3;
rigidificator.ElastificationThreshold = 1e-1; 
integrator = BackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.5, 1e-6 );

energyModel = StVenantKirchoff3DEnergy();

planeContactFinder = PlaneContactFinder3D([0,0,1], [0,0,0], 0.8);
meshMeshContactFinder = MeshSCD3D(0.8,true, []);
contactFinder = {planeContactFinder,meshMeshContactFinder};%

simulate3D(meshes,h,contactFinder, integrator, rigidificator, settings, energyModel);