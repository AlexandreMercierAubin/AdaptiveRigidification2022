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
materials = TriangleMaterial(rho, mu, lambda, alpha0, alpha1);

settings = SimulationSettings();
settings.DrawTimings = 0;
settings.CamPadding(3) = 1;
settings.MakeVideo = 1;
settings.SceneName = 'cantilever_ERPull';
settings.FramesToRecord = 1500;
settings.InitialWindowPosition = [0, 0, 1920, 1080];
% settings.PlotEDotHist = 1;
settings.PlotSkip = 1500;
% settings.WarmStartEnabled = false;
settings.sameComparisonAinv = true;
settings.RecordFramePositionInTD = true;
settings.FirstFrameRigidification = false; % seems to cause a small error for low threshold meshes
% settings.recomputeCacheAinv = true;
scale = [1,1];
rot = 0;

mesh2d = fetchPoly2D('cantileverP05',true, materials,scale, rot, settings);

inds = 1:mesh2d(1).N;
mesh2d.pin(inds(abs(mesh2d.p(1:2:end) - min(mesh2d.p(1:2:end))) < 1e-2));
mesha = AdaptiveMesh(mesh2d);
lenghtMesh = max(mesha.p(1:2:end)) - min(mesha.p(1:2:end));


be = LDLBackwardEuler();
be.Gravity = g;
be.separateQuicksolveGravity = false;

gridX = 10;
gridY = 5;

powerX = logspace(-3,0,gridY);

multY = logspace(-10,-2,gridX);

meshes = {};
rigidificators = {};
testi = 1;
for x = 1:gridY
    for y = 1:gridX
        meshes{testi} = AdaptiveMesh(mesh2d);
%         meshes{testi}.setRigidTransform(0,[5*x,3*y]);%Do not test with
%         this on, otherwise correct the displacement before computing the
%         error
        
        rigidificators{testi} = EDotMexRigidificator();
        rigidificators{testi}.RigidificationThreshold = powerX(x);
        
        %elastification multiplyer
        rigidificators{testi}.ElastificationThreshold = multY(y)*powerX(x);
        rigidificators{testi}.Preconditionner = 'ICHOL';
        rigidificators{testi}.FrameCount = 5;
        rigidificators{testi}.Permutation = 'NONE';
        testi = testi+1;
    end
end
meshes{testi} = mesh2d;
rigidificators{testi} = EDotMexRigidificator();
rigidificators{testi}.RigidificationThreshold = 0;
rigidificators{testi}.ElastificationThreshold = 0;
% rigidificators{testi}.separateQuicksolveGravity = false;

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
save("MEResPull_"+string(gridX)+"_"+string(gridY)+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');

ERres = AnalyseTDER(td, testi, repmat(powerX,gridX,1), repmat(multY,gridY,1)');
threshRatio = ERres(:,3)./ERres(:,4);
color = mat2gray(log10(threshRatio));
figure;
speedFactor = 1./ERres(:,1);
maxError = ERres(:,2);
scatter(speedFactor, log10(maxError), 10, [color,zeros(size(color)),color]);
xlabel("speed factor");
ylabel("log max error scaled by length");
writematrix(ERres,"MEResPull_"+datestr(now,'mm-dd-yyyy_HH-MM')+".txt");
plotTauData(ERres,gridSize);
