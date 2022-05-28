cla;
clear;
close all;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
clear meshes;
rho = 10;
nu = 0.35;      % Poisson ratio: close to 0.5 for rubber
k = 5e4;     % Young's modulus: 0.01e9 approximate for rubber
[mu,lambda] = toLame(nu,k);    % Lamé parameter
alpha0 = 0.00001;   % Rayleigh factor on M
tMaterial1 = [TriangleMaterial(rho, mu, lambda, alpha0, 0.005)];
tMaterial2 = [TriangleMaterial(rho, mu, lambda, alpha0, 0.01)];
tMaterial3 = [TriangleMaterial(rho, mu, lambda, alpha0, 0.05)];
tMaterial4 = [TriangleMaterial(rho, mu, lambda, alpha0, 0.1)];
tMaterial5 = [TriangleMaterial(rho, mu, lambda, alpha0, 0.5)];

s = 6;
dist = 3;

settings.MakeVideo = 1;
settings.FramesToRecord = 2000;
% settings.PGSiterations = 20;
% settings.DrawVelocities = true;
% settings.DrawForces = true;
% settings.DrawLambdas = 1;
% settings.DrawRigidDOFs = 1;
settings.campos=[12,16,8];
settings.projection = "orthographic";
% settings.PlotEDotHist = 1;
settings.PGSiterations = 30;
settings.SceneName = 'beamDamping';
settings.OBJDir = './objs/beamDamping/';
settings.WriteOBJs = true;

resetMeshes = false;
baseMesh = meshLoader("beam913", [], tMaterial5, [0.6,0.6,0.6], resetMeshes, settings);
baseMesh.setRigidTransform([0,90,0],[0,s,1]);
xPos = baseMesh.p(1:3:end);
pinInd = find(xPos == min(xPos));
baseMesh.pin(sort(pinInd));
meshes1 = AdaptiveMesh3D(baseMesh);

baseMesh = meshLoader("beam913", [], tMaterial4, [0.6,0.6,0.6], resetMeshes, settings);
baseMesh.setRigidTransform([0,90,0],[0,s-dist,1]);
baseMesh.pin(sort(pinInd));
meshes2 = AdaptiveMesh3D(baseMesh);

baseMesh = meshLoader("beam913", [], tMaterial3, [0.6,0.6,0.6], resetMeshes, settings);
baseMesh.setRigidTransform([0,90,0],[0,s-2*dist,1]);
baseMesh.pin(sort(pinInd));
meshes3 = AdaptiveMesh3D(baseMesh);

baseMesh = meshLoader("beam913", [], tMaterial2, [0.6,0.6,0.6], resetMeshes, settings);
baseMesh.setRigidTransform([0,90,0],[0,s-3*dist,1]);
baseMesh.pin(sort(pinInd));
meshes4 = AdaptiveMesh3D(baseMesh);

baseMesh = meshLoader("beam913", [], tMaterial1, [0.6,0.6,0.6], resetMeshes, settings);
baseMesh.setRigidTransform([0,90,0],[0,s-4*dist,1]);
baseMesh.pin(sort(pinInd));
meshes5 = AdaptiveMesh3D(baseMesh);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 1e-7;
rigidificator.ElastificationThreshold = 1e-5; 

rigidificator2 = EDot3DMexRigidificator();
rigidificator2.RigidificationThreshold = 1e-20;
rigidificator2.ElastificationThreshold = 1; 

integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.3, 1e-3 );

energyModel = NeoHookean3DEnergy();

nullContactFinder = NullContactFinder(3);
contactFinder = {nullContactFinder};

simulate3D({meshes1,meshes2,meshes3,meshes4,meshes5},h,contactFinder, integrator, {rigidificator,rigidificator,rigidificator,rigidificator,rigidificator}, settings, energyModel);