clear

h = 0.01; % time step

rho = 50;
nu = 0.4; % Poisson ratio: close to 0.5 for rubber
E = 5e3; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.5; % Rayleigh factor on M
alpha1 = 0.0;%1; % Rayleigh factor on K
material = TriangleMaterial(rho, mu, lambda, alpha0, alpha1);

resetMesh = true;
mesh2d =  fetchPoly2D('cantileverp05',resetMesh, material);
%mesh = AdaptiveMesh(loadMeshFromPOLY('data/cantileverp3',material));
mesh2d.setRigidTransform( 90, [0,0] );
% Awkward to pin based on index numbers... but works here!  Could also
% consider pinning DOFs based on if they fall within a region.
mesh2d.pin(mesh2d.N-5:mesh2d.N);

mesh2d = AdaptiveMesh(mesh2d);

settings = SimulationSettings();
% settings.DrawTimings = 0;
settings.RunRightAway = 1;
% settings.DrawForces = 0;
settings.DrawEDots = 1;
settings.CamPadding(2) = 0;
settings.MakeVideo = true;
settings.SceneName = "hangingCantileverHD";
settings.recomputeCacheAinv = true;
% settings.plotTriImplicitRigidificationElastification = 1;
rigi = EDotMexRigidificator();

rigi.FrameCount = 3;
rigi.RigidificationThreshold = 1e-5;
rigi.ElastificationThreshold = 6e-4;
% gshift = GravityChangeScript();d
% gshift.frameNumbers = [400];
% gshift.gravity = [0.05];
integrator = LDLBackwardEuler();
integrator.separateQuicksolveGravity = false;

simulate( mesh2d, integrator, h, settings, rigi, NullContactFinder(), NullContactFinder());%, gshift
