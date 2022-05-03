cla;
clear;
close all;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear meshes;
rho = 40;
nu = 0.22;      % Poisson ratio: close to 0.5 for rubber
k = 5e3;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame(nu, k);
alpha0 = 0.001;   % Rayleigh factor on M
alpha1 = 0.01;  % Rayleigh factor on K
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];
tMaterial.cacheName = 'blobsMat';

resetMeshCache = false;
baseMesh = fetchMesh3D('blobsHD',resetMeshCache, @dittosBuilder, tMaterial);
baseMesh.setRigidTransform([0,0,0],[0,0,0.5]);
extraCantilever = Mesh3D(baseMesh);
extraCantilever.setRigidTransform([0,0,0],[0,0,2]);

baseMesh.mergeMesh(extraCantilever);

meshes = AdaptiveMesh3D(baseMesh);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 1e-5;
rigidificator.ElastificationThreshold = 1e-4; 
rigidificator.Permutation = 'DISSECT';
rigidificator.Preconditionner = 'ICHOLICT';
integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0.2, 0.1 );
integrator.Gravity = -9.8;
% integrator.useFullAinv = true;

energyModel = NeoHookean3DEnergy();
% energyModel = CorotationalEnergy();
% energyModel = StVenantKirchoff3DEnergy();

planeContactFinder = PlaneContactFinder3D([0,0,1], [0,0,0], 0.4);
meshMeshContactFinder = MeshSCD3D(0.2,false, []);
contactFinder = {planeContactFinder,meshMeshContactFinder};

settings.PlotSkip = 2;
settings.SceneName = 'blobs';
% settings.InitialWindowPosition = [0,0,1920,1080];
settings.FramesToRecord = 3500;
% settings.WriteOBJs = true;
settings.OBJDir = './objs/blobs/';
settings.MakeVideo = 1;
settings.campos=[-10,10,3];

td = simulate3D({baseMesh, meshes},h,contactFinder, integrator, rigidificator, settings, energyModel);
save("blob_"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
writeTDcsv(td,"dittos",["_default","_adaptive"]);%baseMesh,
