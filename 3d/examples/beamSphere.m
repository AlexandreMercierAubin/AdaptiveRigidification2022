cla;
clear;
close all;
h = 0.001; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
clear meshes;
rho = 10;
nu = 0.40;      % Poisson ratio: close to 0.5 for rubber
k = 1e4;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu,k );
alpha0 = 0.0001;   % Rayleigh factor on M
alpha1 = 0.001;  % Rayleigh factor on K
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

baseMesh = meshLoader("beam", [], tMaterial, [1,1,1], false, settings);
baseMesh.setRigidTransform([0,0,0],[0,0,1.2]);
extraCantilever = Mesh3D(baseMesh);
extraCantilever.setRigidTransform([0,0,0],[0,0,3.2]);
baseMesh.mergeMesh(extraCantilever);
meshes = AdaptiveMesh3D(baseMesh);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 1e-4;
rigidificator.ElastificationThreshold = 9e-4; 
integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.1, 1e-2 );

energyModel = StVenantKirchoff3DEnergy();

planeContactFinder = PlaneContactFinder3D([0,0,1], [0,0,0], 0.4);
sphereContactFinder = SphereContactFinder(0.5, [5,0,2.8], 0.2);
meshMeshContactFinder = MeshSCD3D(0.5,false, []);
contactFinder = {planeContactFinder,sphereContactFinder,meshMeshContactFinder};

settings.SceneName = 'beamSphere';
settings.MakeVideo = 1;
settings.WriteOBJs = 0;
settings.FramesToRecord = 800;

settings.campos=[-17,20,5];

simulate3D(meshes,h,contactFinder, integrator, rigidificator, settings, energyModel);