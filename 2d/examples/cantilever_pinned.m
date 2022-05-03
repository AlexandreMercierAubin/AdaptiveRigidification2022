clear();
close all;
g = -9.8; % gravity
h = 0.01; % time step

rho = 25;
nu = 0.4; % Poissonff ratio: close to 0.5 for rubber
E = 1e5; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.01; % Rayleigh factor on M
alpha1 = 0.01; % Rayleigh factor on K
material = TriangleMaterial(rho, mu, lambda, alpha0, alpha1);

settings = SimulationSettings();
settings.PlotEDotHist = 1;
settings.SceneName = "pinnedCantileverLowDef";
settings.MakeVideo = 0;
% settings.FramesToRecord = 700;
settings.DrawDv = true;
settings.DrawApproxDv = true;
settings.PCGiterations = 1; % under 5 the dv start spiraling
% settings.quicksolveSimulation = true;

rot = 0;
scale = [1,1];
mesh2d = fetchPoly2D('cantileverP3',true, material, scale, rot, settings);

xPos = mesh2d.p(1:2:end);
pinInd = find(xPos <= min(xPos)+0.01);
mesh2d.pin(sort(pinInd));

be = LDLBackwardEuler();
be.Gravity = g;
be.separateQuicksolveGravity = false;

rigid = EDotMexRigidificator();
rigid.RigidificationThreshold = 1e-6;
rigid.ElastificationThreshold = 5e-6;

mesha = AdaptiveMesh(mesh2d);

simulate( mesha, be, h, settings, rigid );
