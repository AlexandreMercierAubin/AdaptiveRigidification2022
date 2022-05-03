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

baseMesh = meshLoader("beam730", [], tMaterial, [0.6,0.6,0.6], false, settings);
baseMesh.setRigidTransform([0,90,0],[0,0,0]);
zPos = baseMesh.p(3:3:end);
bottomInds = find(zPos <= min(zPos) + 0.05);
dofs = reshape([bottomInds*3],1,[])';

yPos = baseMesh.p(2:3:end);
farRightInds = find(yPos <= min(yPos) + 0.05);
dofs2 = reshape([farRightInds*3-1],1,[])';

animScripted = SequentialPositionAnimationScripter();
animScripted.frameNumbers = [1,2,3,70]; %force position of the bottom nodes on those frames
animScripted.dofs = {dofs,dofs,dofs,dofs2};
p1 = baseMesh.p(dofs) +0.05;
p2 = baseMesh.p(dofs) +0.1;
p3 = baseMesh.p(dofs) +0.15;
p4 = baseMesh.p(dofs2) +0.5;
animScripted.positions = {p1, p2, p3, p4};

meshes = AdaptiveMesh3D(baseMesh);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 1e-3;
rigidificator.ElastificationThreshold = 1e-1; 
integrator = BackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.2, 1e-4 );

energyModel = StVenantKirchoff3DEnergy();

nullContactFinder = NullContactFinder(3);
contactFinder = {nullContactFinder};

settings.MakeVideo = 1;
settings.FramesToRecord = 200;
% settings.PGSiterations = 20;
% settings.DrawVelocities = true;
% settings.DrawForces = true;
settings.DrawLambdas = 1;
settings.campos=[10,10,1];
settings.PGSiterations = 30;
settings.SceneName = 'beamJump';
simulate3D(meshes,h,contactFinder, integrator, rigidificator, settings, energyModel, animScripted);