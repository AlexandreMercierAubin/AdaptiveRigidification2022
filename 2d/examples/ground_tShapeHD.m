clear;

h = 0.01; % time step

rho = 10;
nu = 0.4; % Poisson ratio: close to 0.5 for rubber
E = 1e5; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.03;
alpha1 = 0.001;

material1 = TriangleMaterial(rho, mu, lambda, alpha0, alpha1);

material2 = TriangleMaterial(rho*10, mu, lambda, alpha0, alpha1);
material2.color = [ .4 .8 .3 ];

materials = [ material1, material2 ];



m1 = fetchPoly2D('TShapeHD.1',true, materials);

% compute the element centers to adjust masses on part of the T
v = reshape(m1.p,2,[]);
n1 = v( :, m1.t(:,1) );
n2 = v( :, m1.t(:,2) );
n3 = v( :, m1.t(:,3) );
elCenter = 1/3 * (n1+n2+n3);
% if y > .4 set to heavy material
ind = elCenter(2,:) > 0.4;
attributes = m1.materialIndex;
attributes(ind) = 2;
m1.updateMaterials( attributes, materials );

m2 = AdaptiveMesh( m1 );

m1.setRigidTransform( 0, [-.75 0] );
m2.setRigidTransform( 0, [.75 0] );

settings = SimulationSettings();
settings.DrawTimings = 0;
settings.CamPadding(3) = 3;
settings.DrawEdges = true;
settings.SceneName = "ground_tShapeHD";
settings.DrawLambdasScale = 1;
settings.DrawLambdas = 1;
settings.WarmStartEnabled = true;
settings.PGSiterations = 300;

rigid = EDotMexRigidificator(); 
rigid.RigidificationThreshold = 1e-5;
rigid.ElastificationThreshold = 1e-4;

integrator = LDLBackwardEuler();
integrator.Gravity = -9.8; 
integrator.separateQuicksolveGravity = true;
comp = 1e-3;
ERP = 0.1;%.1;
% integrator.setComplianceAndBaumgarteFromERPandCFM(h,1, 0.001 );
% integrator.fullElasticContacts = 1;
integrator.Baumgarte = ERP/h;
integrator.Compliance = comp ;
pcf = PlaneContactFinder([0.0, 1.0] , [0, -1], 0);

% comparison!
meshes = {m1, m2};
td = simulate( meshes, integrator, h, settings, rigid, pcf );


