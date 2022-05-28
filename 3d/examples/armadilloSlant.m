cla;
clear;
close all;
h = 0.0001; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear meshes;

rho = 3;
nu = 0.45;      % Poisson ratio: close to 0.5 for rubber
k = 1e3;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame(nu, k);
alpha0 = 0.000000;   % Rayleigh factor on M
alpha1 = 0.00000;  % Rayleigh factor on K
tMaterial1 = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

rho = 3;
nu = 0.4;     % Poisson ratio: close to 0.5 for rubber
E = 1e9;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, E );
alpha0 = 0.0000001;%.0001;   % Rayleigh factor on M
alpha1 = 0.1;  % Rayleigh factor on K
tMaterial2 = TriangleMaterial(rho, mu, lambda, alpha0, alpha1,[0.7,0.5,0.7]);

tMaterial = [tMaterial1];
resetMeshCache = false;

settings.MakeVideo = 1;
settings.SceneName = 'armadilloSlant';
settings.FramesToRecord = 50000;
settings.PlotSkip = 100;
% settings.DrawRigidDOFs = true;
% settings.WriteOBJs = true;
settings.OBJDir = './objs/armadilloSlant/';
settings.campos=[7,10,1];
settings.FocusOnMeshNode = [1,30];

baseMesh = fetchMesh3D('armadilloSlant',resetMeshCache, @armadilloSlantBuilder, [tMaterial], settings);
mesha = AdaptiveMesh3D(baseMesh);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 1e-6;
rigidificator.ElastificationThreshold = 5e-6; 
integrator = LDLBackwardEuler3D();
integrator.Gravity = -9.8;
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0.1, 0.1 );

energyModel = NeoHookean3DEnergy();

planeContactFinder = PlaneContactFinder3D([0.2,0.1,0.70], [0,0,0], 0.72);%[0.35,0.15,0.5]
meshMeshContactFinder = MeshSCD3D(0.2,false, []);
contactFinder = {planeContactFinder,meshMeshContactFinder};

simulate3D(mesha,h,contactFinder, integrator, rigidificator, settings, energyModel);