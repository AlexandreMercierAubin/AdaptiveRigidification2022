clear 

h = 0.01; % time step

rho = 10;
nu = 0.30; % Poisson ratio: close to 0.5 for rubber
E = 1e5; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.00;%3;
alpha1 = 0.00;%2;

settings = SimulationSettings();
% settings.MaximizeWindow = 1;
settings.DrawTimings = 0;
%settings.DrawRigidFrames = 0;
%settings.DrawContact = 1;
settings.DrawLambdas = 0;
settings.PrintTimings = 0; % Turning this on is genearlly a bad idea
settings.PGSiterations = 30;

% settings.RunRightAway = false;
settings.CamPadding = [ 5 5 5 5 ];
settings.CamPadding = [ 1 1 1 1 ];

settings.SceneName = "bagelRotate";

material = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1)];
scale = [1,1];
rot = 0;
mesh2d = fetchPoly2D('bagel.1',false, material,scale,rot,settings);

m1a = AdaptiveMesh( mesh2d );


rigid = EDiffMexRigidificator();
rigid.RigidificationThreshold = 1e-5;
rigid.ElastificationThreshold = 1e-4; 
rigid.FrameCount = 3; %number of frames before rigidifying.

integrator = LDLBackwardEuler();
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0.1, 0.001 );

integrator.Gravity = -9.8; % gravity does not seem to be the problem


pcf = PlaneContactFinder([0.0, 1.0] , [0.01, -10], 1);
%ccf = MovingCircleContactFinder( 5.1 , [0, 0], 4 );     % just a bit bigger than the hole
ccf = MovingCircleContactFinder( 0.7 , [0, 0], 4 );     % just a bit smaller than the hole

            ccf.cfun = @(t) [ 0, 0.5*(cos(t*8)-1) * (mod( t*8, pi*2*10 ) < pi*2) ];
            ccf.dcdt = @(t) [ 0, -0.5*sin(t*8)*8  * (mod( t*8, pi*2*10 ) < pi*2) ];

             ccf.thetafun = @(t)  1.5*sin(t)  * (mod( t, pi*2*3 ) < pi*2);
             ccf.dthetadt = @(t)  1.5*cos(t)  * (mod( t, pi*2*3 ) < pi*2);
            
%            ccf.thetafun = @(t)  0.4*t;
%            ccf.dthetadt = @(t)  0.4;

ncf = NullContactFinder;
%cfs = { pcf, ccf };
%cfs = { ncf };
cfs = { ccf };

%meshes = m1a;
td = simulate( m1a, integrator, h, settings, rigid, cfs );
