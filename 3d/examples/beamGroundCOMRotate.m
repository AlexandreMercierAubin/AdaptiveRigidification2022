cla;
clear;
close all;
h = 0.005; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear meshes;
rho = 10;
nu = 0.35;      % Poisson ratio: close to 0.5 for rubber
k = 1e4;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, k );
alpha0 = 0.0001;   % Rayleigh factor on M
alpha1 = 0.001;  % Rayleigh factor on K
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

% settings.MakeVideo = 1;
% settings.FramesToRecord = 600;
settings.SceneName = 'beamCOMRotate';
% settings.PlotEDotHist = 1;
settings.campos=[12,16,8];
settings.projection = "orthographic";
% settings.DrawRigidDOFs = 1;

scale = [1,1,1];
resetMesh = false;

baseMesh = meshLoader("beam", [], tMaterial, scale, resetMesh, settings);
baseMesh.setRigidTransform([0,90,0],[0,-2,1],true);
meshes = AdaptiveMesh3D(baseMesh);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 5e-4;
rigidificator.ElastificationThreshold = 9e-4; 
integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0.1,0.001 );

energyModel = StVenantKirchoff3DEnergy();

planeContactFinder = PlaneContactFinder3D([0,0,1], [0,0,-3], 0.2);
contactFinder = {planeContactFinder};



simulate3D(meshes,h,contactFinder, integrator, rigidificator, settings, energyModel);