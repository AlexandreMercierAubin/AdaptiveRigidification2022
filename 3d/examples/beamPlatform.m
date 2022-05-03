cla;
clear;
close all;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();
% settings.MakeVideo = 1;
% settings.FramesToRecord = 600;
settings.SceneName = 'beamGround';
% settings.PlotEDotHist = 1;
settings.campos=[-17,20,5];
% settings.DrawRigidDOFs = 1;

% MESHES

clear meshes;
rho = 10;
nu = 0.48;      % Poisson ratio: close to 0.5 for rubber
k = 1e4;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, k );
alpha0 = 0.0001;   % Rayleigh factor on M
alpha1 = 0.001;  % Rayleigh factor on K
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

resetMesh = false;
baseMesh = meshLoader("beam", [], tMaterial, [1,1,1], resetMesh, settings);
R = 
meshes = AdaptiveMesh3D(baseMesh);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 5e-4;
rigidificator.ElastificationThreshold = 9e-4; 
integrator = BackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.3, 1e-3 );

energyModel = StVenantKirchoff3DEnergy();


platformContactDetector = PlatformContactDetector([0,0,0], [15,15,1], [0,0,-3], 0.1);
contactFinder = {platformContactDetector};



simulate3D(meshes,h,contactFinder, integrator, rigidificator, settings, energyModel);