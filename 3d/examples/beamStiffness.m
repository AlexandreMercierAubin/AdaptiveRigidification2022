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
k = 1e3;     % Young's modulus: 0.01e9 approximate for rubber
[mu1,lambda1] = toLame(nu,1e4);    % Lamé parameter
[mu2,lambda2] = toLame(nu,3e4);    % Lamé parameter
[mu3,lambda3] = toLame(nu,5e4);    % Lamé parameter
[mu4,lambda4] = toLame(nu,1e5);    % Lamé parameter
[mu5,lambda5] = toLame(nu,1e6);    % Lamé parameter
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.1;  % Rayleigh factor on K
tMaterial1 = [TriangleMaterial(rho, mu1, lambda1, alpha0, alpha1)];
tMaterial2 = [TriangleMaterial(rho, mu2, lambda2, alpha0, alpha1)];
tMaterial3 = [TriangleMaterial(rho, mu3, lambda3, alpha0, alpha1)];
tMaterial4 = [TriangleMaterial(rho, mu4, lambda4, alpha0, alpha1)];
tMaterial5 = [TriangleMaterial(rho, mu5, lambda5, alpha0, alpha1)];

s = 6;
dist = 3;

baseMesh = meshLoader("beam913", [], tMaterial5, [0.6,0.6,0.6], false, settings);
baseMesh.setRigidTransform([0,90,0],[0,s,1]);
xPos = baseMesh.p(1:3:end);
pinInd = find(xPos == min(xPos));
baseMesh.pin(sort(pinInd));
meshes1 = AdaptiveMesh3D(baseMesh);

baseMesh = meshLoader("beam913", [], tMaterial4, [0.6,0.6,0.6]);
baseMesh.setRigidTransform([0,90,0],[0,s-dist,1]);
baseMesh.pin(sort(pinInd));
meshes2 = AdaptiveMesh3D(baseMesh);

baseMesh = meshLoader("beam913", [], tMaterial3, [0.6,0.6,0.6]);
baseMesh.setRigidTransform([0,90,0],[0,s-2*dist,1]);
baseMesh.pin(sort(pinInd));
meshes3 = AdaptiveMesh3D(baseMesh);

baseMesh = meshLoader("beam913", [], tMaterial2, [0.6,0.6,0.6]);
baseMesh.setRigidTransform([0,90,0],[0,s-3*dist,1]);
baseMesh.pin(sort(pinInd));
meshes4 = AdaptiveMesh3D(baseMesh);

baseMesh = meshLoader("beam913", [], tMaterial1, [0.6,0.6,0.6]);
baseMesh.setRigidTransform([0,90,0],[0,s-4*dist,1]);
baseMesh.pin(sort(pinInd));
meshes5 = AdaptiveMesh3D(baseMesh);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 1e-7;
rigidificator.ElastificationThreshold = 1e-5; 

rigidificatorSoft = EDot3DMexRigidificator();
rigidificatorSoft.RigidificationThreshold = 1e-6;
rigidificatorSoft.ElastificationThreshold = 1e-4; 

integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.3, 1e-3 );

energyModel = NeoHookean3DEnergy();

nullContactFinder = NullContactFinder(3);
contactFinder = {nullContactFinder};

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
settings.SceneName = 'beamStiffness';
settings.OBJDir = './objs/beamStiffness/';
settings.WriteOBJs = true;
meshes = {meshes1,meshes2,meshes3,meshes4,meshes5};
rigidificators = {rigidificatorSoft,rigidificatorSoft,rigidificator,rigidificator,rigidificator};
simulate3D(meshes,h,contactFinder, integrator, rigidificators, settings, energyModel);