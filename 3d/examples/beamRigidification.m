cla;
clear;
close all;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();
%NOT A VERY USEFUL TEST TBH
% MESHES
clear meshes;
rho = 10;
nu = 0.35;      % Poisson ratio: close to 0.5 for rubber
k = 1e5;     % Young's modulus: 0.01e9 approximate for rubber
[mu1,lambda1] = toLame(nu,k);    % Lamé parameter
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.001;  % Rayleigh factor on K
tMaterial1 = [TriangleMaterial(rho, mu1, lambda1, alpha0, alpha1)];

baseMesh = meshLoader("beam730", [], tMaterial1, [0.6,0.6,0.6], false, settings);
baseMesh.setRigidTransform([0,90,0],[0,1,1]);
xPos = baseMesh.p(1:3:end);
pinInd = find(xPos == min(xPos));
baseMesh.pin(sort(pinInd));
meshes = AdaptiveMesh3D(baseMesh);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 1e-7;
rigidificator.ElastificationThreshold = 1; 

rigidificator2 = EDot3DMexRigidificator();
rigidificator2.RigidificationThreshold = 1e-6;
rigidificator2.ElastificationThreshold = 1; 

rigidificator3 = EDot3DMexRigidificator();
rigidificator3.RigidificationThreshold = 1e-5;
rigidificator3.ElastificationThreshold = 1; 

rigidificator4 = EDot3DMexRigidificator();
rigidificator4.RigidificationThreshold = 1e-4;
rigidificator4.ElastificationThreshold = 1; 

integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.3, 1e-3 );

energyModel = NeoHookean3DEnergy();

planeContactFinder = NullContactFinder(3);
contactFinder = {planeContactFinder};

settings.MakeVideo = 1;
settings.FramesToRecord = 3000;
% settings.PGSiterations = 20;
% settings.DrawVelocities = true;
% settings.DrawForces = true;
% settings.DrawLambdas = 1;
% settings.DrawRigidDOFs = 1;
settings.campos=[20,0,1];
% settings.PlotEDotHist = 1;
settings.PGSiterations = 30;
settings.SceneName = 'beamRigidification';
settings.OBJDir = './objs/beamRigidification/';
settings.WriteOBJs = true;
simulate3D(meshes,h,contactFinder, integrator, {rigidificator,rigidificator2,rigidificator3,rigidificator4}, settings, energyModel);