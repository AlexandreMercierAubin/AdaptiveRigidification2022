cla;
clear;
close all;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear meshes;

rho = 5;
nu = 0.35;      % Poisson ratio: close to 0.5 for rubber
k = 5e3;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, k );
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.1;  % Rayleigh factor on K
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [1,1,1])];

rho = 7;
nu = 0.35;      % Poisson ratio: close to 0.5 for rubber
k = 5e4;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, k );
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.05;  % Rayleigh factor on K
tMaterial2 = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [1,0,0])];

materials = [tMaterial,tMaterial2];

baseMesh = fetchMesh3D('pachinkoPill',false, @pachinkoPillHeterogenousBuilder, materials, settings);
meshes = AdaptiveMesh3D(baseMesh);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 5e-3;
rigidificator.ElastificationThreshold = 5e-2; 
rigidificator.FrameCount = 7;
integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0.2,0.01 );
%apparently the ERP needs to be super big in this scene otherwise objects
%fall through the pachinko
% energyModel = StVenantKirchoff3DEnergy();
energyModel = NeoHookean3DEnergy();

rotation1 = [0,-25,0];
rotation2 = [0,25,0];
scale = [15, 15, 1];

platformContactDetector1 = ObjectContactDetector("3d/data/pachinko1Model.obj",[90,0,0], [20,12,6], [-3,3,0],0.05);
meshMeshContactFinder = MeshSCD3D(0.1, true, []);
contactFinder = {platformContactDetector1,meshMeshContactFinder};

settings.MakeVideo = 1;
settings.FramesToRecord = 2600;
settings.SceneName = 'PachinkoPillHeterogeneous';
% settings.PlotEDotHist = 1;
settings.PlotSkip = 1;
settings.campos=[300,0,5];
settings.camtarget = [0,0,3];
% settings.FocusOnMeshNode = [1,30];
settings.camfov = 10;
settings.WriteOBJs = 1;
% settings.DrawRigidDOFs = 1;
settings.OBJDir = './objs/pachinkoPill/';

td = simulate3D({baseMesh,meshes},h,contactFinder, integrator, rigidificator, settings, energyModel);%baseMesh,
save("pills_"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
writeTDcsv(td,"PachinkoPill",["_default", "_adaptive"]);%"_default",