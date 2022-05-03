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
k = 1e5;     % Young's modulus: 0.01e9 approximate for rubber
[mu1,lambda1] = toLame(0.30,k);    % Lamé parameter
[mu2,lambda2] = toLame(0.35,k);    % Lamé parameter
[mu3,lambda3] = toLame(0.4,k);    % Lamé parameter
[mu4,lambda4] = toLame(0.45,k);    % Lamé parameter
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.1;  % Rayleigh factor on K
tMaterial1 = [TriangleMaterial(rho, mu1, lambda1, alpha0, alpha1)];
tMaterial2 = [TriangleMaterial(rho, mu2, lambda2, alpha0, alpha1)];
tMaterial3 = [TriangleMaterial(rho, mu3, lambda3, alpha0, alpha1)];
tMaterial4 = [TriangleMaterial(rho, mu4, lambda4, alpha0, alpha1)];

baseMesh = meshLoader("beam730", [], tMaterial4, [0.6,0.6,0.6], false, settings);
baseMesh.setRigidTransform([0,90,0],[0,4,1]);
xPos = baseMesh.p(1:3:end);
pinInd = find(xPos == min(xPos));
baseMesh.pin(sort(pinInd));
meshes4 = AdaptiveMesh3D(baseMesh);

baseMesh = meshLoader("beam730", [], tMaterial3, [0.6,0.6,0.6], false, settings);
baseMesh.setRigidTransform([0,90,0],[0,2.5,1]);
baseMesh.pin(sort(pinInd));
meshes3 = AdaptiveMesh3D(baseMesh);

baseMesh = meshLoader("beam730", [], tMaterial2, [0.6,0.6,0.6], false, settings);
baseMesh.setRigidTransform([0,90,0],[0,1,1]);
baseMesh.pin(sort(pinInd));
meshes2 = AdaptiveMesh3D(baseMesh);

baseMesh = meshLoader("beam730", [], tMaterial1, [0.6,0.6,0.6], false, settings);
baseMesh.setRigidTransform([0,90,0],[0,-0.5,1]);
baseMesh.pin(sort(pinInd));
meshes1 = AdaptiveMesh3D(baseMesh);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 1e-7;
rigidificator.ElastificationThreshold = 1e-5; 
integrator = BackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.3, 1e-3 );

energyModel = NeoHookean3DEnergy();

nullContactFinder = NullContactFinder(3);
contactFinder = {nullContactFinder};

settings.MakeVideo = 1;
settings.FramesToRecord = 3000;
% settings.PGSiterations = 20;
% settings.DrawVelocities = true;
% settings.DrawForces = true;
% settings.DrawLambdas = 1;
% settings.DrawRigidDOFs = 1;
settings.campos=[12,16,8];
settings.projection = "orthographic";
% settings.PlotEDotHist = 1;
settings.PGSiterations = 30;
settings.SceneName = 'beamPoisson';
settings.OBJDir = './objs/beamPoisson/';
settings.WriteOBJs = true;
simulate3D({meshes1,meshes2,meshes3,meshes4},h,contactFinder, integrator, rigidificator, settings, energyModel);