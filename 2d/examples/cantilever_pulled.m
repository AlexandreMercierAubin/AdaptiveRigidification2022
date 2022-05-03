clear();
close all;
g = -9.8; % gravity
h = 0.01; % time step

rho = 30;
nu = 0.4; % Poisson ratio: close to 0.5 for rubber
E = 1e5; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.0001; % Rayleigh factor on M
alpha1 = 0.05; % Rayleigh factor on K
material = TriangleMaterial(rho, mu, lambda, alpha0, alpha1);

settings = SimulationSettings();
settings.DrawTimings = 0;
settings.CamPadding(3) = 3;
% settings.MakeVideo = 1;
settings.SceneName = 'cantilever_ERPull';
settings.FramesToRecord = 1500;
% settings.PlotEDotHist = 1;
% settings.PlotSkip = 1500;
% settings.WarmStartEnabled = false;
settings.sameComparisonAinv = true;
settings.RecordFramePositionInTD = true;
% settings.recomputeCacheAinv = true;

scale = [1,1];
rot = 0;

mesh2d = fetchPoly2D('cantileverP05',true, material,scale, rot, settings);

inds = 1:mesh2d(1).N;
mesh2d.pin(inds(abs(mesh2d.p(1:2:end) - min(mesh2d.p(1:2:end))) < 1e-2));
mesha = AdaptiveMesh(mesh2d);
lenghtMesh = max(mesha.p(1:2:end)) - min(mesha.p(1:2:end));


be = LDLBackwardEuler();
be.Gravity = g;
be.separateQuicksolveGravity = false;

gridSize = 7;
powerX = 1e10^(1/(gridSize-1));

multY = logspace(1,10,gridSize);

meshes = {};
rigidificators = {};

meshes{1} = AdaptiveMesh(mesh2d);
rigidificators = EDotMexRigidificator();
rigidificators.FrameCount = 5;
rigidificators.RigidificationThreshold = 1e-4;
rigidificators.ElastificationThreshold = 1e-3;
rigidificators.Permutation = 'NONE';
rigidificators.Preconditionner = 'ICHOL';

meshes{2} = mesh2d;

n1 = 700;
n2 = 725;
vertInds = inds(abs(mesh2d.p(1:2:end) - max(mesh2d.p(1:2:end))) < 1e-2);
dofs = vertInds*2;

animScript = ForceImpulseAnimationScripter();
animScript.frameNumbers = n1:n2-1; %force position of the bottom nodes on those frames
animScript.dofs = repelem({dofs},n2-n1);
initForce = zeros(numel(dofs),1);
animScript.forceImpulse = repelem({initForce},n2-n1);
%Drops by 0.1 for n2-n1 frames
forceImpulse = 2;
for i = 2:n2-n1
    animScript.forceImpulse{i}(1:end) = animScript.forceImpulse{i-1}(1:end) + forceImpulse;
end

td = simulate( meshes, be, h, settings, rigidificators, NullContactFinder, NullContactFinder, animScript);

