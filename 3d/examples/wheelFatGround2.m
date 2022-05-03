% cla;
% clear;
close all;
h = 0.001; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear meshes;

%tire
rho = 3;
nu = 0.45;      % Poisson ratio: close to 0.5 for rubber
k = 5e3;   % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame(nu, k);
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.05;  % Rayleigh factor on K
tMaterial1 = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

%rim
rho = 10;
nu = 0.4;     % Poisson ratio: close to 0.5 for rubber
E = 1e6;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, E );
alpha0 = 0.00001;%.0001;   % Rayleigh factor on M
alpha1 = 0.01;  % Rayleigh factor on K
tMaterial2 = TriangleMaterial(rho, mu, lambda, alpha0, alpha1,[0.7,0.5,0.7]);

tMaterial = [tMaterial1,tMaterial2];

% baseMesh = eleNodeLoader("wheelFat-regions.1", [], tMaterial, [1,1,1]);
baseMesh = meshLoader("wheelFatHD",  [], tMaterial, [1,1,1], false, settings);
v = reshape(baseMesh.p,3,[]);
n1 = v( :, baseMesh.t(:,1) );
n2 = v( :, baseMesh.t(:,2) );
n3 = v( :, baseMesh.t(:,3) );
n4 = v( :, baseMesh.t(:,4) );
elCenter = 0.25 * (n1+n2+n3+n4);

dfromcenter = sqrt(sum(elCenter.^2, 1));
inds = dfromcenter < 0.5;
attributes = baseMesh.materialIndex;
attributes(inds) = 2;
baseMesh.updateMaterials( attributes, tMaterial );

baseMesh.setRigidTransform([270,0,90],[0,0,1.1]);
baseMesh.setRigidMotion([0,2,0],[1,0,0]);
mesh = Mesh3D(baseMesh);
mesha = AdaptiveMesh3D(baseMesh);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 5e-4;
rigidificator.ElastificationThreshold = 5e-3; 
rigidificator.Preconditionner = 'ICHOLICT';
rigidificator.Permutation = 'DISSECT';
integrator = LDLBackwardEuler3D();
integrator.Gravity = -9.8;
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0.2, 0.1);

energyModel = NeoHookean3DEnergy();

planeContactFinder = PlaneContactFinder3D([0,0,1], [0,0,0], 0.72);
contactFinder = {planeContactFinder};

settings.MakeVideo = 1;
settings.SceneName = 'wheelFatGround2';
settings.FramesToRecord = 7000;
% settings.PlotSkip = 10000;
% settings.DrawRigidDOFs = true;
% settings.WriteOBJs = true;
settings.OBJDir = './objs/wheelFatGround2/';
settings.campos=[4,7,1]*1.5;
settings.camtarget = [0.5,0,1.5];
settings.FocusOnMeshNode = [1,30];
% settings.PlotEDotHist = 1;
td = simulate3D({mesh,mesha},h,contactFinder, integrator, rigidificator, settings, energyModel);
save("wheel_"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
% writeTDcsv(td,"wheelFatGround2",["_default","_adaptive"]);%mesh,