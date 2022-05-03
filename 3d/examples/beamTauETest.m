clear();
close all;
g = -9.8; % gravity
h = 0.01; % time step

rho = 10;
nu = 0.35;      % Poisson ratio: close to 0.5 for rubber
k = 5e4;     % Young's modulus: 0.01e9 approximate for rubber
[mu,lambda] = toLame(nu,k);    % Lamé parameter
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.1;  % Rayleigh factor on K
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];
tMaterial.cacheName = 'ScalingTestMat';


settings = Simulation3DSettings();
settings.DrawTimings = 0;
settings.CamPadding(3) = 3;
% settings.MakeVideo = 1;
settings.SceneName = 'beamTauETest';
settings.FramesToRecord = 500;
% settings.PlotEDotHist = 1;
% settings.PlotSkip = 500;
% settings.WarmStartEnabled = false;
% settings.PGSiterations = 300;
settings.sameComparisonAinv = true;
settings.RecordFramePositionInTD = true;
settings.recomputeCacheAinv = true;
settings.campos=[20,0,1];

resetMesh = false;
scale = [1,1];
rot = 0;
baseMesh5232 = meshLoader("beam5232", [], tMaterial, [0.6,0.6,0.6], false, settings);
% baseMesh.setRigidTransform([0,90,0],[0,s,1]);
baseMesh5232.setRigidTransform([0,90,0],[0,0,0]);
xPos = baseMesh5232.p(1:3:end);
pinInd = find(xPos == min(xPos));
baseMesh5232.pin(sort(pinInd));

be = LDLBackwardEuler3D();
be.Gravity = g;
be.separateQuicksolveGravity = false;
be.setComplianceAndBaumgarteFromERPandCFM(h,0.2,0.1);

energyModel = NeoHookean3DEnergy();
nullContactFinder = NullContactFinder(3);
contactFinder = {nullContactFinder};

gridSize = 20;

powerX = logspace(-10,0,gridSize);

meshes = {};
rigidificators = {};
testi = 1;
for x = 1:gridSize
    meshes{testi} = AdaptiveMesh3D(baseMesh5232);
%         meshes{testi}.setRigidTransform(0,[5*x,3*y]);%Do not test with
%         this on, otherwise correct the displacement before computing the
%         error

    rigidificators{testi} = EDot3DMexRigidificator();
    rigidificators{testi}.RigidificationThreshold = powerX(x)*0.1;

    %elastification multiplyer
    rigidificators{testi}.ElastificationThreshold = powerX(x);
    rigidificators{testi}.Preconditionner = 'ICHOLICT';
    rigidificators{testi}.FrameCount = 3;
    rigidificators{testi}.Permutation = 'DISSECT';
    testi = testi+1;
end
meshes{testi} = baseMesh5232;
rigidificators{testi} = EDot3DMexRigidificator();
rigidificators{testi}.RigidificationThreshold = 0;
rigidificators{testi}.ElastificationThreshold = 0;
% rigidificators{testi}.separateQuicksolveGravity = false;

td = simulate3D( meshes, h, nullContactFinder, be, rigidificators, settings, energyModel);
save("beamTauERes_"+string(gridSize)+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
AnalyseTDtauE(td, powerX, 3);

