cla;
clear;
close all;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear meshes;
%for octopusP2 resolution
rho = 10;
nu = 0.4;      % Poisson ratio: close to 0.5 for rubber
k = 3e3;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, k );
alpha0 = 0.0001;   % Rayleigh factor on M
alpha1 = 0.1;  % Rayleigh factor on K
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];

%for octopusP1
% rho = 3;
% nu = 0.4;      % Poisson ratio: close to 0.5 for rubber
% k = 5e3;     % Young's modulus: 0.01e9 approximate for rubber
% [ mu, lambda ] = toLame( nu, k );
% alpha0 = 0.0001;   % Rayleigh factor on M
% alpha1 = 0.08;  % Rayleigh factor on K
% tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];


resetMeshCache = false;
baseMesh = fetchMesh3D('octopus',resetMeshCache, @octopusBowlBuilder, [tMaterial]);
meshes = AdaptiveMesh3D(baseMesh);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 3e-5;
rigidificator.ElastificationThreshold = 3e-4; 
rigidificator.Preconditionner = 'ICHOLICT';
rigidificator.Permutation = "DISSECT";
integrator = LDLBackwardEuler3D();
integrator.separateQuicksolveGravity = true;
% integrator.PGSquicksolve = false;
% integrator.useFullAinv = true;
% integrator.regularizator = 5;
settings.recomputeCacheAinv = true;
comp = 0.01;
ERP = 0.2;

integrator.setComplianceAndBaumgarteFromERPandCFM(h,ERP,comp);

% energyModel = StVenantKirchoff3DEnergy();
energyModel = NeoHookean3DEnergy();

friction = 0.5;
objContactDetector = ObjectContactDetector("3d/data/fatBowl.obj",[90,0,0], [5,5,5], [0,0,0],friction);
contactFinder = {objContactDetector};

% settings.MakeVideo = 1;
settings.FramesToRecord = 4000;
settings.SceneName = 'octopusBowl';
% settings.WriteOBJs = 1;
settings.OBJDir = './objs/octopusBowlICT/';
settings.campos=[15,15,5];
% settings.PlotSkip = 1;
% settings.PlotSkip = 4000;
settings.camtarget = [0,0,0];
settings.PGSiterations = 600;
% settings.recomputeCacheAinv = true;
% settings.PCGiterations = 1;
% settings.quicksolveSimulation = false;


animScripted = ForceImpulseAnimationScripter();
start = 1100;
entries = 40;
xPos = baseMesh.p(1:3:end);
tuggedInds = find(xPos >= max(xPos) - 0.8);
dofs = reshape([tuggedInds*3-2,tuggedInds*3-1,tuggedInds*3]',[],1)';
 
impulse = 55;
forceFloat = zeros(size(dofs));
forceFloat(3:3:end) = meshes.mass(dofs(3:3:end))*(-integrator.Gravity+impulse);
forceFloat(2:3:end) = -meshes.mass(dofs(2:3:end))*impulse;
forceFloat(1:3:end) = meshes.mass(dofs(1:3:end))*impulse;

for i = 1:entries
    animScripted.dofs{i} = dofs;
    animScripted.forceImpulse{i} = forceFloat;
end
animScripted.frameNumbers = start:start+entries-1;
%baseMesh, 
td = simulate3D({baseMesh,meshes},h,contactFinder, integrator, rigidificator, settings, energyModel,animScripted);
save("octopusPlatform_"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
writeTDcsv(td,"octopusPlatform",["_default","_adaptive"]);%baseMesh,