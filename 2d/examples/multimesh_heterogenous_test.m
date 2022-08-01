clear all;
close all;

g = -9.8; % gravity
h = 0.01; % time step

ERP = 0.1;
CFM = 0.01;

friction = 0.2;

rho = 5;
nu1 = 0.4; % Poisson ratio: close to 0.5 for rubber
E1 = 1e4; % Young's modulus: 0.01e9 approximate for rubber
[ mu1, lambda1 ] = toLame( nu1, E1 );

rho2 = 1;
nu2 = 0.4; % Poisson ratio: close to 0.5 for rubber
E2 = 5e3; % Young's modulus: 0.01e9 approximate for rubber
[ mu2, lambda2 ] = toLame( nu2, E2 );

alpha01 = 0.01; % Rayleigh factor on M
alpha11 = 0.05;    % Rayleigh factor on K
alpha02 = 0.01; % Rayleigh factor on M
alpha12 = 0.05; % Rayleigh factor on K
materials = [TriangleMaterial(rho, mu1, lambda1, alpha01, alpha12, [0.1,0.7,0.1]),... 
             TriangleMaterial(rho2, mu2, lambda2, alpha02, alpha12, [0.1,0.1,0.7])]

m1 =  fetchPoly2D('cantileverP3',true, materials);
v = reshape(m1.p,2,[]);
n1 = v( :, m1.t(:,1) );
n2 = v( :, m1.t(:,2) );
n3 = v( :, m1.t(:,3) );
elCenter = (n1+n2+n3)/3;
ends = find(elCenter(1,:) <= -0.5 | elCenter(1,:) >= 0.5);
attributes = m1.materialIndex;
attributes(ends) = 2;
m1.updateMaterials( attributes, materials );

m2 = Mesh(m1);
m3 = Mesh(m1);

m1.setRigidTransform( 0, [0,1.6] );
m2.setRigidTransform( 0, [0,3] );
m2.setRigidMotion( 0, [-2.5 -1] );
m3.setRigidTransform( 0, [0,4.2] );

mesh = m1;
mesh.mergeMesh( m2 );
mesh.mergeMesh( m3 );
mesha = AdaptiveMesh( mesh );

settings = SimulationSettings();
% settings.MakeVideo = 1;
% settings.FramesToRecord = 500;
%settings.MaximizeWindow = 1;
settings.DrawTimings = 0;
settings.CamPadding(3) = 5;
% settings.DrawContact = 1;
settings.DrawLambdas = 1;
settings.ElasticContacts = 1;
settings.PGSiterations = 50;
settings.recomputeCacheAinv = true;
settings.WarmStartEnabled = false;
% settings.PCGiterations = 500;
% settings.quicksolveSimulation = 1;

settings.SceneName = 'multimerge_heterogenous';
rigid = EDotMexRigidificator();
% rigid = EDotVectorRigidificator();
rigid.FrameCount = 3;
rigid.RigidificationThreshold = 1e-4;
rigid.ElastificationThreshold = 1e-2;

integrator = BackwardEuler();
integrator.regularizator = 6;

integrator.setComplianceAndBaumgarteFromERPandCFM(h,ERP, CFM );
pcf = PlaneContactFinder([0.0, 1.0] , [0, -1], friction );
mcf = MeshSCD( friction, false );
td = simulate( mesha, integrator, h, settings, rigid, pcf, mcf );
