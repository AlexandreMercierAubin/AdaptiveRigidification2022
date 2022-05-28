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
[mu1,lambda1] = toLame(nu,k);    % Lamé parameter
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.1;  % Rayleigh factor on K
tMaterial1 = [TriangleMaterial(rho, mu1, lambda1, alpha0, alpha1)];

dist = 3;
resetMeshes = false;

% settings.MakeVideo = 1;
settings.FramesToRecord = 2000;
% settings.PGSiterations = 20;
% settings.DrawVelocities = true;
% settings.DrawForces = true;
% settings.DrawLambdas = 1;
% settings.DrawRigidDOFs = 1;
settings.campos=[12,16,8];
settings.projection = "orthographic";
settings.PGSiterations = 30;
settings.SceneName = 'beamElastification';
settings.OBJDir = './objs/beamElastification/';
% settings.WriteOBJs = true;

baseMesh = meshLoader("beam913", [], tMaterial1, [0.6,0.6,0.6], resetMeshes, settings);
baseMesh.setRigidTransform([0,90,0],[0,dist,1]);
xPos = baseMesh.p(1:3:end);
pinInd = find(xPos == min(xPos));
baseMesh.pin(sort(pinInd));
meshes1 = AdaptiveMesh3D(baseMesh);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 1e-8;
rigidificator.ElastificationThreshold = 1e-6; 

baseMesh.setRigidTransform([0,0,0],[0,-dist,0]);
meshes2 = AdaptiveMesh3D(baseMesh);
rigidificator2 = EDot3DMexRigidificator();
rigidificator2.RigidificationThreshold = 1e-7;
rigidificator2.ElastificationThreshold = 1e-5; 

baseMesh.setRigidTransform([0,0,0],[0,-dist,0]);
meshes3 = AdaptiveMesh3D(baseMesh);
rigidificator3 = EDot3DMexRigidificator();
rigidificator3.RigidificationThreshold = 1e-6;
rigidificator3.ElastificationThreshold = 1e-4; 

baseMesh.setRigidTransform([0,0,0],[0,-dist,0]);
meshes4 = AdaptiveMesh3D(baseMesh);
rigidificator4 = EDot3DMexRigidificator();
rigidificator4.RigidificationThreshold = 1e-5;
rigidificator4.ElastificationThreshold = 1e-3; 

notUsed = EDot3DMexRigidificator();
notUsed.RigidificationThreshold = 1e-4;
notUsed.ElastificationThreshold = 1e-3; 

baseMesh.setRigidTransform([0,0,0],[0,dist*4,0]);

integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.3, 1e-3 );

energyModel = NeoHookean3DEnergy();

planeContactFinder = NullContactFinder(3);
contactFinder = {planeContactFinder};

rigidificators ={rigidificator,rigidificator2,rigidificator3,rigidificator4, notUsed};
meshes = {meshes1,meshes2,meshes3,meshes4, baseMesh};
simulate3D(meshes,h,contactFinder, integrator, rigidificators, settings, energyModel);
l2m1 = sqrt(sum((baseMesh.p - meshes1.p).^2))/baseMesh.N
l2m2 = sqrt(sum((baseMesh.p - meshes2.p).^2))/baseMesh.N
l2m3 = sqrt(sum((baseMesh.p - meshes3.p).^2))/baseMesh.N
l2m4 = sqrt(sum((baseMesh.p - meshes4.p).^2))/baseMesh.N
pb = reshape(baseMesh.p, 3, size(baseMesh.p, 1) / 3)';

p1 = reshape(meshes1.p, 3, size(meshes1.p, 1) / 3)';
p1(:,2) = p1(:,2) +3;
p2 = reshape(meshes2.p, 3, size(meshes2.p, 1) / 3)';
p2(:,2) = p2(:,2) +6;
p3 = reshape(meshes3.p, 3, size(meshes3.p, 1) / 3)';
p3(:,2) = p3(:,2) +9;
p4 = reshape(meshes4.p, 3, size(meshes4.p, 1) / 3)';
p4(:,2) = p4(:,2) +12;
diff1 = max(sqrt(sum((pb-p1).^2,2)))
diff2 = max(sqrt(sum((pb-p2).^2,2)))
diff3 = max(sqrt(sum((pb-p3).^2,2)))
diff4 = max(sqrt(sum((pb-p4).^2,2)))
