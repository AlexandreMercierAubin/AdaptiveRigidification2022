%cla;
clear;
%close all;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear mesh3d;
rho = 10;
nu = 0.35;      % Poisson ratio: close to 0.5 for rubber
k = 4e4;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, k );
alpha0 = 0.0001;   % Rayleigh factor on M
alpha1 = 0.001;  % Rayleigh factor on K
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];
tMaterial.cacheName = 'beamGroundMat';

% % settings.MakeVideo = 1;
%settings.FramesToRecord = 1000;
settings.SceneName = 'beamGround';
% settings.PlotEDotHist = 1;
settings.campos=[12,16,8]*1.3;
%settings.projection = "orthographic";
% settings.DrawRigidDOFs = 1;

baseMesh = meshLoader("beam", [], tMaterial,[1,1,1],false, settings);
normal = [0.2,0,0.8];
R = findVectorRotation(normal);
baseMesh.setRigidRotationFromMatrix(R,true);
mesh3d = AdaptiveMesh3D(baseMesh);
mesh3d.setRigidTransform([0,0,0],[0,0,3]);
% meshesichol = AdaptiveMesh3D(baseMesh);
% meshesichol.setRigidTransform([0,0,0],[0,3,0]);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 5e-5;
rigidificator.ElastificationThreshold = 1e-3;  

integrator = LDLBackwardEuler3D();
ERP = 0.2;
CFM = 0.001;
integrator.setComplianceAndBaumgarteFromERPandCFM( h, ERP, CFM );

rigidificatoriChol = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 1e-5;
rigidificator.ElastificationThreshold = 1e-3; 
rigidificator.Preconditionner = 'ICHOLICT';

energyModel = NeoHookean3DEnergy();

position = [0 0 -1];
frictionCoefficient = 0.5;
planeContactFinder = PlaneContactFinder3D( normal, position, frictionCoefficient );
contactFinder = {planeContactFinder};

% m = {mesh3d,meshesichol};
% r = {rigidificator,rigidificatoriChol}

m = {mesh3d};
r = {rigidificator};
simulate3D( m, h, contactFinder, integrator, r, settings, energyModel );

