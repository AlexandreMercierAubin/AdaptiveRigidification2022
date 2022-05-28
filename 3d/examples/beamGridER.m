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

% settings.MakeVideo = 1;
settings.FramesToRecord = 2000;
settings.campos=[12,16,8];
settings.projection = "orthographic";
settings.PGSiterations = 30;
settings.SceneName = 'beamElastification';
settings.OBJDir = './objs/beamElastification/';
settings.PlotSkip = 2000;
% settings.WriteOBJs = true;
settings.sameComparisonAinv = true;

baseMesh = meshLoader("beam913", [], tMaterial1, [0.6,0.6,0.6], false, settings);
lenghtMesh = max(baseMesh.p(3:3:end)) - min(baseMesh.p(3:3:end));
baseMesh.setRigidTransform([0,90,0],[0,0,1]);
xPos = baseMesh.p(1:3:end);
pinInd = find(xPos == min(xPos));
baseMesh.pin(sort(pinInd));
meshes1 = AdaptiveMesh3D(baseMesh);

gridSize = 1;
meshes = {};
rigidificators = {};
testi = 1;
for i = 1:gridSize*gridSize
    x = floor((i-1)/gridSize);
    y = mod((i-1),gridSize);
    if x > y
        continue;
    end
    meshes{i} = AdaptiveMesh3D(baseMesh);
    rigidificators{i} = EDot3DMexRigidificator();
    RT = 1e-10*power^x;
    rigidificators{i}.RigidificationThreshold = RT;
    ET = 1e-10*power^y;
    rigidificators{i}.ElastificationThreshold = ET;
    testi = testi+1;
end
meshes{testi} = baseMesh;

%non adaptive mesh so whatever
rigidificators{testi} = EDot3DMexRigidificator();
rigidificators{testi}.RigidificationThreshold = 0;
rigidificators{testi}.ElastificationThreshold = 0; 

integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.1, 1e-3 );

energyModel = NeoHookean3DEnergy();

planeContactFinder = NullContactFinder(3);
contactFinder = {planeContactFinder};

td = simulate3D(meshes,h,contactFinder, integrator, rigidificators, settings, energyModel);
pb = reshape(meshes{testi}.p, 3, size(meshes{testi}.p, 1) / 3)';

testLast = testi-1;
maxErr = zeros(testLast,1);
speedup = zeros(testLast,1);
timeEl = sum(td{testLast}.log(2,:));
tre = zeros(testLast,1);
trr = zeros(testLast,1);
for i = 1:testLast
    pi = reshape(meshes{i}.p, 3, size(meshes{i}.p, 1) / 3)';
    maxErr(i) = max(sqrt(sum((pb-pi).^2,2)))/lenghtMesh;
    time = sum(td{i}.log(2,:));
    speedup(i) = timeEl/time;
    tre(i) = rigidificators{i}.ElastificationThreshold;
    trr(i) = rigidificators{i}.RigidificationThreshold;
end

ERresGrid3D = [speedup,maxErr,tre,trr];
writematrix(ERresGrid3D);

tauE = log10(ERresGrid3D(:,3));
tauR = log10(ERresGrid3D(:,4));
su = ERresGrid3D(:,1);
err = ERresGrid3D(:,2);

set(gcf,'color','w');

ix = tauE >= tauR;

subplot(1,2,1); 
scatter(tauR(ix),tauE(ix),5,(err(ix)),'filled');
axis equal
colorbar
xlim([-10,0]);
ylim([-10,0]);
xlabel('tauR');
ylabel('tauE');
title( sprintf("max error\n") );

subplot(1,2,2); 
scatter(tauR(ix),tauE(ix),5,su(ix),'filled');
axis equal
colorbar
xlim([-10,0]);
ylim([-10,0]);
xlabel('tauR');
ylabel('tauE');
title(sprintf("speedup\n") );
