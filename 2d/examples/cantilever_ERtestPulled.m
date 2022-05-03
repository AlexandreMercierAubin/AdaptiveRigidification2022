clear();
close all;
g = -9.8; % gravity
h = 0.01; % time step

rho = 30;
nu = 0.4; % Poisson ratio: close to 0.5 for rubber
E = 1e5; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.0001; % Rayleigh factor on M
alpha1 = 0.01; % Rayleigh factor on K
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

gridSize = 10;
power = 1e10^(1/(gridSize-1));
meshes = {};
rigidificators = {};
testi = 1;
for i = 1:gridSize*gridSize
    x = floor((i-1)/gridSize);
    y = mod((i-1),gridSize);
    if x > y
        continue;
    end
    meshes{testi} = AdaptiveMesh(mesh2d);
    rigidificators{testi} = EDotMexRigidificator();
    RT = 1e-10*power^x;
    rigidificators{testi}.RigidificationThreshold = RT;
    ET = 1e-10*power^y;
    rigidificators{testi}.ElastificationThreshold = ET;
    testi = testi+1;
end
meshes{testi} = mesh2d;
rigidificators{testi} = EDotMexRigidificator();
rigidificators{testi}.RigidificationThreshold = 0;
rigidificators{testi}.ElastificationThreshold = 0;

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
for i = 2:n2-n1
    animScript.forceImpulse{i}(1:end) = animScript.forceImpulse{i-1}(1:end) + 1;
end

td = simulate( meshes, be, h, settings, rigidificators, NullContactFinder, NullContactFinder, animScript);

for frame = 1:size(td{testi}.p,2)
        pb{frame} = reshape(td{testi}.p(:,frame), 2, size(meshes{testi}.p, 1) / 2)';
end

maxErr = zeros(testi-1,1);
sumErr = zeros(testi-1,1);
speedup = zeros(testi-1,1);
timeEl = sum(td{testi}.log(2,:));
tre = zeros(testi-1,1);
trr = zeros(testi-1,1);
for i = 1:testi-1
    for frame = 1:size(td{i}.p,2)
        pi = reshape(td{i}.p(:,frame), 2, size(meshes{i}.p, 1) / 2)';
        errDistance = sum(sqrt(sum((pb{frame}-pi).^2,2)));
        sumErr(i) = sumErr(i)+ errDistance;
        maxErrTmp = max(errDistance)/lenghtMesh;
        maxErr(i) = max(maxErr(i),maxErrTmp);
    end
    time = sum(td{i}.log(2,:));
    speedup(i) = timeEl/time;
    tre(i) = rigidificators{i}.ElastificationThreshold;
    trr(i) = rigidificators{i}.RigidificationThreshold;
end
ERres = [speedup,maxErr,tre,trr];
writematrix(ERres,"ERes_"+"ERes_"+datestr(now,'mm-dd-yyyy_HH-MM')+".txt"+".txt");

tauE = log10(ERres(:,3));
tauR = log10(ERres(:,4));
su = ERres(:,1);
err = log10(ERres(:,2));

set(gcf,'color','w');

ix = tauE>=tauR;

subplot(1,2,1); 
scatter(tauR(ix),tauE(ix),5,(err(ix)),'filled');
axis equal
c = colorbar;
xlim([-10,0]);
ylim([-10,0]);
xlabel('tauR');
ylabel('tauE');
title( sprintf("max error\n") );
c.Label.String = 'log 10 scale';

subplot(1,2,2); 
scatter(tauR(ix),tauE(ix),5,su(ix),'filled');
axis equal
c = colorbar;
xlim([-10,0]);
ylim([-10,0]);
xlabel('tauR');
ylabel('tauE');
title(sprintf("speedup\n") );
c.Label.String = 'log 10 scale';

%fig 2
figure;
ERres = [speedup,sumErr,tre,trr];
writematrix(ERres);

tauE = log10(ERres(:,3));
tauR = log10(ERres(:,4));
su = ERres(:,1);
err = log10(ERres(:,2));

set(gcf,'color','w');

ix = tauE>=tauR;

subplot(1,2,1); 
scatter(tauR(ix),tauE(ix),5,(err(ix)),'filled');
axis equal
c = colorbar;
xlim([-10,0]);
ylim([-10,0]);
xlabel('tauR');
ylabel('tauE');
title( sprintf("sum error distance\n") );
c.Label.String = 'log 10 scale';

subplot(1,2,2); 
scatter(tauR(ix),tauE(ix),5,su(ix),'filled');
axis equal
c = colorbar;
xlim([-10,0]);
ylim([-10,0]);
xlabel('tauR');
ylabel('tauE');
title(sprintf("speedup\n") );
c.Label.String = 'log 10 scale';

figure;
scatter(tauR(ix),tauE(ix),5,ERres(ix,2)./su(ix),'filled');
c = colorbar;
xlim([-10,0]);
ylim([-10,0]);
xlabel('tauR');
ylabel('tauE');
title( sprintf("error/speedup\n") );

figure;
scatter(tauR(ix),tauE(ix),5,err(ix)./su(ix),'filled');
c = colorbar;
xlim([-10,0]);
ylim([-10,0]);
xlabel('tauR');
ylabel('tauE');
title( sprintf("log10 error/speedup\n") );