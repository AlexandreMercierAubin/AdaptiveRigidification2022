clear();
g = -9.8; % gravity
h = 0.01; % time step

rho = 30;
nu = 0.4; % Poisson ratio: close to 0.5 for rubber
E = 1e5; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0; % Rayleigh factor on M
alpha1 = 0.01; % Rayleigh factor on K

% Set damping from damping factor
% compute first modal frequency via eigs
% where omega0 is sqrt(D(1,1))
% (set a breakpoint in LDLBackwardEuler to harvest this)
% Mii=M(ii,ii);Kii=(K(ii,ii)+K(ii,ii)')/2;
% [V,D] = eigs(Kii,Mii,7,'smallestabs'); 
% omeag0 = sqrt(-D(1,1))
% If unpinned 2D then use D(4,4)
% If unpinned 3D then use D(7,7)
omega0 = 6.6; % divide by 2pi to get period
df = 0.2;  % zero is undamped, 1 is critical
% set alpha1 from alpha0
alpha0 = 0;
alpha1 = (2*df-alpha0/omega0)/omega0;
% set alpha0 from alpha1
% alpha1 = 0;
% alpha0 = (2*df-alpha1*omega0)*omega0;

materials = TriangleMaterial(rho, mu, lambda, alpha0, alpha1);

settings = SimulationSettings();
settings.DrawTimings = 0;
settings.CamPadding(3) = 3;
settings.MakeVideo = 0;
settings.SceneName = 'cantilever_mega';
settings.FramesToRecord = 1000;
% settings.PlotEDotHist = 1;
% settings.WarmStartEnabled = false;
% settings.PGSiterations = 100;
settings.plotTriImplicitRigidificationElastification = 1;


resetMesh = false;
scale = [1,1];
rot = 0;
mesh2d = fetchPoly2D('cantileverP05',resetMesh, materials, scale, rot, settings);

inds = 1:mesh2d(1).N;
mesh2d.pin(inds(abs(mesh2d.p(1:2:end) - min(mesh2d.p(1:2:end))) < 1e-2));
mesha = AdaptiveMesh(mesh2d);

rigid = EDotMexRigidificator();
% rigid = EDotVectorRigidificator();
rigid.RigidificationThreshold = 1e-6;
rigid.ElastificationThreshold = 1e-0;
% rigid.ScaleByMaxEdgeLength = 1;

be = LDLBackwardEuler();
be.Gravity = g;

td = simulate( {mesh2d,mesha}, be, h, settings, rigid );

if ( numel(td) > 1 ) 
    figure(4); clf;
    td{1}.plotCompareSiggraphVersion( td{2} ,h);
    pbaspect([3 2 1])
    title('Simulation Breakdown', 'FontName', 'Linux Biolinum O');
    xlabel('Time(s)', 'FontName', 'Linux Biolinum O');
    ylabel('Time(ms)', 'FontName', 'Linux Biolinum O');
    saveas(gcf,'out/simulation.pdf');
    
    figure(6); clf;
    plot(mesha.N-td{2}.logCounts(2,4:end),'Color', [204/255, 85/255, 0], 'LineWidth', 2);
    title('Number of Rigid Particles per Frame', 'FontName', 'Linux Biolinum O');
    xlabel('Time(s)', 'FontName', 'Linux Biolinum O');
    ylabel('number of particles', 'FontName', 'Linux Biolinum O');
    pbaspect([3 2 1]);
    saveas(gcf,'out/rigidParticles.pdf');

    figure(7); clf;
    plot(td{2}.logCounts(7,4:end),'Color', [204/255, 85/255, 0], 'LineWidth', 2 );
    title('Number of Degrees of Freedom per Frame', 'FontName', 'Linux Biolinum O');
    xlabel('Time(s)', 'FontName', 'Linux Biolinum O');
    ylabel('DOFs', 'FontName', 'Linux Biolinum O');
    pbaspect([3 2 1]);
    saveas(gcf,'out/totalDofs.pdf');
end

