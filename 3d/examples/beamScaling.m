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
[mu,lambda] = toLame(nu,k);    % Lamé parameter
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.1;  % Rayleigh factor on K
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];
tMaterial.cacheName = 'ScalingTestMat';

s=6;
dist = 3;

% baseMesh28 = meshLoader("beam28", [], tMaterial, [0.6,0.6,0.6], false, settings);
% % baseMesh.setRigidTransform([0,90,0],[0,s-dist*2,1]);
% baseMesh28.setRigidTransform([0,90,0],[0,0,0]);
% xPos = baseMesh28.p(1:3:end);
% pinInd = find(xPos == min(xPos));
% baseMesh28.pin(sort(pinInd));
% meshes28 = AdaptiveMesh3D(baseMesh28);
% 
% baseMesh365 = meshLoader("beam365", [], tMaterial, [0.6,0.6,0.6], false, settings);
% % baseMesh.setRigidTransform([0,90,0],[0,s-dist,1]);
% baseMesh365.setRigidTransform([0,90,0],[0,0,0]);
% xPos = baseMesh365.p(1:3:end);
% pinInd = find(xPos == min(xPos));
% baseMesh365.pin(sort(pinInd));
% meshes365 = AdaptiveMesh3D(baseMesh365);
% 
% baseMesh913 = meshLoader("beam913", [], tMaterial, [0.6,0.6,0.6], false, settings);
% % baseMesh.setRigidTransform([0,90,0],[0,s-dist,1]);
% baseMesh913.setRigidTransform([0,90,0],[0,0,0]);
% xPos = baseMesh913.p(1:3:end);
% pinInd = find(xPos == min(xPos));
% baseMesh913.pin(sort(pinInd));
% meshes913 = AdaptiveMesh3D(baseMesh913);
% 
baseMesh1264 = meshLoader("beam1264", [], tMaterial, [0.6,0.6,0.6], false, settings);
% baseMesh.setRigidTransform([0,90,0],[0,s-dist,1]);
baseMesh1264.setRigidTransform([0,90,0],[0,0,0]);
xPos = baseMesh1264.p(1:3:end);
pinInd = find(xPos == min(xPos));
baseMesh1264.pin(sort(pinInd));
meshes1264 = AdaptiveMesh3D(baseMesh1264);

baseMesh2248 = meshLoader("beam2248", [], tMaterial, [0.6,0.6,0.6], false, settings);
% baseMesh.setRigidTransform([0,90,0],[0,s-dist,1]);
baseMesh2248.setRigidTransform([0,90,0],[0,0,0]);
xPos = baseMesh2248.p(1:3:end);
pinInd = find(xPos == min(xPos));
baseMesh2248.pin(sort(pinInd));
meshes2248 = AdaptiveMesh3D(baseMesh2248);

baseMesh5232 = meshLoader("beam5232", [], tMaterial, [0.6,0.6,0.6], false, settings);
% baseMesh.setRigidTransform([0,90,0],[0,s,1]);
baseMesh5232.setRigidTransform([0,90,0],[0,0,0]);
xPos = baseMesh5232.p(1:3:end);
pinInd = find(xPos == min(xPos));
baseMesh5232.pin(sort(pinInd));
meshes5232 = AdaptiveMesh3D(baseMesh5232);

% baseMesh8325 = meshLoader("beam8325", [], tMaterial, [0.6,0.6,0.6], false, settings);
% % baseMesh.setRigidTransform([0,90,0],[0,s,1]);
% baseMesh8325.setRigidTransform([0,90,0],[0,0,0]);
% xPos = baseMesh8325.p(1:3:end);
% pinInd = find(xPos == min(xPos));
% baseMesh8325.pin(sort(pinInd));
% meshes8325 = AdaptiveMesh3D(baseMesh8325);
% 
% baseMesh16106 = meshLoader("beam16106", [], tMaterial, [0.6,0.6,0.6], false, settings);
% % baseMesh.setRigidTransform([0,90,0],[0,s,1]);
% baseMesh16106.setRigidTransform([0,90,0],[0,0,0]);
% xPos = baseMesh16106.p(1:3:end);
% pinInd = find(xPos == min(xPos));
% baseMesh16106.pin(sort(pinInd));
% meshes16106 = AdaptiveMesh3D(baseMesh16106);
% 
% baseMesh44608 = meshLoader("beam44608", [], tMaterial, [0.6,0.6,0.6], false, settings);
% % baseMesh.setRigidTransform([0,90,0],[0,s,1]);
% baseMesh44608.setRigidTransform([0,90,0],[0,0,0]);
% xPos = baseMesh16106.p(1:3:end);
% pinInd = find(xPos == min(xPos));
% baseMesh44608.pin(sort(pinInd));
% meshes44608 = AdaptiveMesh3D(baseMesh44608);

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 1e-7;
rigidificator.ElastificationThreshold = 3e-1; 
integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.3, 1e-3 );

energyModel = NeoHookean3DEnergy();

nullContactFinder = NullContactFinder(3);
contactFinder = {nullContactFinder};

% settings.MakeVideo = 1;
settings.FramesToRecord = 500;
% settings.PGSiterations = 20;
% settings.DrawVelocities = true;
% settings.DrawForces = true;
% settings.DrawLambdas = 1;
% settings.DrawRigidDOFs = 1;
settings.campos=[12,16,8];
settings.projection = "orthographic";
% settings.PlotEDotHist = 1;
settings.PGSiterations = 30;
settings.SceneName = 'beamScaling';
settings.OBJDir = './objs/beamScaling/';
% settings.WriteOBJs = true;
settings.recomputeCacheAinv = true;

td = simulate3D({baseMesh1264,meshes1264,baseMesh2248,meshes2248,baseMesh5232,meshes5232},h,contactFinder, integrator, rigidificator, settings, energyModel);%baseMesh28,meshes28,baseMesh365,meshes365,baseMesh913,meshes913,baseMesh1264,meshes1264,baseMesh2248,meshes2248,baseMesh5232,meshes5232,baseMesh8325,meshes8325,baseMesh16106,meshes16106,baseMesh44608,meshes44608
writeTDcsv(td, "scaling", ["_el1264","_1264","_el2248","_2248","_el5232","_5232"]);%["_el28","_28","_el365","_365","_el913","_913","_el1264","_1264","_el2248","_2248","_el5232","_5232","_el8325","_8325","_el16106","_16106","_el44608","_44608"]
time = compareSpeedFactors(td)