cla; 
clear;
close all;
h = 0.005; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear meshes;

%pinetree
rho = 50;
nu = 0.37;      % Poisson ratio: close to 0.5 for rubber
k = 2e5;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, k );
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.01;  % Rayleigh factor on K
woodcolor = [145, 132, 240]./255;
tMaterialPine = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1,woodcolor)];
tMaterialPine(1).cacheName = 'PineTree';

%tree
rho = 100;
nu = 0.37;      % Poisson ratio: close to 0.5 for rubber
k = 7e5;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, k );
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.07;  % Rayleigh factor on K
woodcolor = [145, 132, 240]./255;
tMaterialTree = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1,woodcolor)];
tMaterialTree(1).cacheName = 'Tree';

%leaves
rho = 0.1;
nu = 0.37;      % Poisson ratio: close to 0.5 for rubber
k = 7e5;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, k );
alpha0 = 0.0001;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
leavescolor = [132, 240, 199]./255;
tMaterialLeaves = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1,leavescolor)];
tMaterialLeaves(1).cacheName = 'Leaves';

%axe
rho = 45;
nu = 0.37;      % Poisson ratio: close to 0.5 for rubber
k = 5e4;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, k );
alpha0 = 0.0001;   % Rayleigh factor on M
alpha1 = 0.1;  % Rayleigh factor on K
Baseballcolor= [227, 240, 132]./255;
tMaterialAxe = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, Baseballcolor)];
tMaterialAxe(1).cacheName = 'Axe';


baseMesh = fetchMesh3D('forest',false, @forestBuilder, [tMaterialTree,tMaterialPine,tMaterialAxe,tMaterialLeaves], settings);
axeRanges = load('3d/data/cached/forestAxeRange');
axeRange = axeRanges.axeRanges{1};
axeTriRange = axeRanges.axeRanges{2};

meshes = AdaptiveMesh3D(baseMesh);
meshes.ForceElastifyElementInds = [axeTriRange];

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 5e-5;
rigidificator.ElastificationThreshold = 5e-2; 
rigidificator.Preconditionner = 'ICHOLICT';
integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0.2, 0.01 );
integrator.separateQuicksolveGravity = false;

energyModel = StVenantKirchoff3DEnergy();

meshContactFinder = MeshSCD3D(0.8,true, []);
contactFinder = {meshContactFinder};

% settings.MakeVideo = 1;
settings.FramesToRecord = 1000;
settings.SceneName = 'forest';
% settings.PlotEDotHist = 1;
settings.campos=[-4,15,6];
settings.camtarget = [0,0,2];
% settings.DrawRigidDOFs = 1;
settings.OBJDir = './objs/forest/';
% settings.WriteOBJs = true;
% settings.PlotSkip = 1000;

animScripted = ForceImpulseAnimationScripter();

dofs = axeRange'; 
% verts = find(meshes.p(x) > max(meshes.p(x))-1);
% dofs = dofs(reshape([verts*3-2,verts*3-1, verts*3]',1,[])');

posInit = meshes.p(dofs);
entries = 100;% hang in the air
forceFloat = zeros(size(dofs));
forceFloat(3:3:end) = -meshes.mass(dofs(3:3:end))*integrator.Gravity;
forceFloatAmortized = zeros(size(dofs));
forceFloatAmortized(3:3:end) = -meshes.mass(dofs(3:3:end))*(integrator.Gravity+2);

for i = 1:entries
    animScripted.dofs{i} = dofs;
    animScripted.forceImpulse{i} = forceFloat;
end

lastEntries = entries;
entries = lastEntries + 50;
impulse = 30;
leftPush = zeros(size(dofs));
leftPush(1:3:end) = meshes.mass(dofs(1:3:end))*impulse;
for i = lastEntries+1:entries
    animScripted.dofs{i} = dofs;
    animScripted.forceImpulse{i} = forceFloat + leftPush;
end

lastEntries = entries;
entries = lastEntries + 100;
rightPush = zeros(size(dofs));
impulse = 30;
rightPush(1:3:end) = -meshes.mass(dofs(1:3:end))*impulse;
for i = lastEntries+1:entries
    animScripted.dofs{i} = dofs;
    animScripted.forceImpulse{i} = forceFloatAmortized + rightPush;
end

lastEntries = entries;
entries = lastEntries + 10;

for i = lastEntries+1:entries
    animScripted.dofs{i} = dofs;
    animScripted.forceImpulse{i} = forceFloat + leftPush*0.25;
end


lastEntries = entries;
entries = lastEntries + 150;
forwardPush = zeros(size(dofs));
impulse = 20;
forwardPush(2:3:end) = -meshes.mass(dofs(2:3:end))*impulse;
for i = lastEntries+1:entries
    animScripted.dofs{i} = dofs;
    animScripted.forceImpulse{i} = forceFloat + forwardPush + leftPush*0.4;
end

lastEntries = entries;
entries = lastEntries + 100;
for i = lastEntries+1:entries
    animScripted.dofs{i} = dofs;
    animScripted.forceImpulse{i} = forceFloat + leftPush*0.25;
end

animScripted.frameNumbers = 1:entries;

% start = entries+1;
% endEntries = start + 20;
% for i = start:endEntries
%     animScripted.dofs{i} = dofs;
%     posInit(1:3:end) = posInit(1:3:end) - 0.0001;
%     animScripted.positions{i} = posInit;
% end
% 
% start = endEntries+1;
% endEntries = start + 20;
% for i = start:endEntries
%     animScripted.dofs{i} = dofs;
%     posInit(1:3:end) = posInit(1:3:end) + 0.0001;
%     animScripted.positions{i} = posInit;
% end


td = simulate3D({baseMesh,meshes},h,contactFinder, integrator, rigidificator, settings, energyModel, {animScripted});
save("forest_"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
writeTDcsv(td, "forest", ["_default","_adaptive"]);
