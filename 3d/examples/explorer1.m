clear;

h = 0.001; % time step

settings = Simulation3DSettings();

% MESHES

%body
rho = 10;
nu = 0.35;      % Poisson ratio: close to 0.5 for rubber
E = 5e6;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, E );
alpha0 = 0.000001;   % Rayleigh factor on M

alpha1 = 0.00001;  % Rayleigh factor on K
tMaterial1 = TriangleMaterial(rho, mu, lambda, alpha0, alpha1,[0.7,0.5,0.7]);

%anthena
rho = 10;
nu = 0.35;      % Poisson ratio: close to 0.5 for rubber
E = 2e5;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, E );
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.01;%001;  % Rayleigh factor on K
tMaterial2 = TriangleMaterial(rho, mu, lambda, alpha0, alpha1,[0.5,0.7,0.7]);

%end of the antena
rho = 15;
nu = 0.35;      % Poisson ratio: close to 0.5 for rubber
E = 1e7;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, E );
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.001;%001;  % Rayleigh factor on K
tMaterial3 = TriangleMaterial(rho, mu, lambda, alpha0, alpha1,[0.7,0.7,0.2]);

%soft planet
rho = 15;
nu = 0.35;      % Poisson ratio: close to 0.5 for rubber
E = 1e4;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame( nu, E );
alpha0 = 0.0001;   % Rayleigh factor on M
alpha1 = 0.001;%001;  % Rayleigh factor on K
tMaterialplanet = TriangleMaterial(rho, mu, lambda, alpha0, alpha1,[0.7,0.7,0.7]);

materials = [ tMaterial1, tMaterial2, tMaterial3];

mesh = meshLoader("explorer1", [], materials, [0.1,0.1,0.1], false, settings);
mesh.setRigidTransform( [90,0,0],[0,0,0.3]);
mesh.setRigidMotion([0,15,0],[0,0,0]);

mesh2 = meshLoader("PlanetHigher2", [], tMaterialplanet, [1,1,1], false, settings);
dist = [4,-13,3];
mesh2.setRigidTransform( [0,0,0],dist*3);
mesh2.setRigidMotion([0,0,0.5],-dist);

% mesh3 = meshLoader("THandle", [], materials, [0.1,0.1,0.1]);
% mesh3.setRigidTransform( [0,0,0], [0,3,3]); % put it higher in the world frame
% mesh3.setRigidMotion( [0,-5,0], [0,0,0] );

% some care to make this a multi-material model...
% let's first figure out where the centers of the elements are
v = reshape(mesh.p,3,[]);
n1 = v( :, mesh.t(:,1) );
n2 = v( :, mesh.t(:,2) );
n3 = v( :, mesh.t(:,3) );
n4 = v( :, mesh.t(:,4) );
elCenter = 0.25 * (n1+n2+n3+n4);
% The .24 here is because it doesn't quite sit on the y axis... small shift
% in the z direction.
dfromyaxis = sqrt(sum((elCenter([1,3],:)-[0;.24]).^2, 1));
rthresh = 0.7;
% find the bits away from the core, and avoid the boxes on the side of teh
% rocket or other features that are too far off the centerline
ind = (dfromyaxis > rthresh) & (abs(elCenter(2,:))<1);
% Use a scatter to make sure we got the rigth bits!
%scatter3( elCenter(1,ind), elCenter(2,ind), elCenter(3,ind),'.')

attributes = mesh.materialIndex;
attributes(ind) = 2;

ind = (dfromyaxis > 3.1) & (abs(elCenter(2,:))<1);
attributes(ind) = 3;

mesh.updateMaterials( attributes, materials );

% Give the antennas a tug
dvertfromyaxis = sqrt(sum((v([1,3],:)-[0;.24]).^2, 1));
vind = (dvertfromyaxis > rthresh) & (abs(v(2,:))<1);
vy = mesh.v(2:3:end);
onx = abs(v(1,:))<1;
vy(vind) = vy(vind) + 2*(dvertfromyaxis(vind)' - rthresh) .* (1 - 2*onx(vind)');
mesh.v(2:3:end) = vy;

mesh.mergeMesh(mesh2);
% mesh.mergeMesh(mesh3);
mesha = AdaptiveMesh3D(mesh);

% pinning tris
% zPos = meshes.p(3:3:end);
% min_I = find(zPos(:,1) <= min(zPos)+2);
% meshes.pin(sort(min_I));

rigidificator = EDot3DMexRigidificator();
rigidificator.RigidificationThreshold = 1e-4;
rigidificator.ElastificationThreshold = 1e-3; 
rigidificator.FrameCount = 3;
rigidificator.Preconditionner = 'ICHOLICT';
rigidificator.Permutation = 'DISSECT';

integrator = LDLBackwardEuler3D();
integrator.Gravity = 0;
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0.2, 0.01 );
energyModel = NeoHookean3DEnergy();

meshMeshContactFinder = MeshSCD3D(0.5,true, []);
contactFinder = {meshMeshContactFinder};

% CONFIG

settings.RunRightAway = 1;

settings.PlotSkip = 10;

settings.campos = [-300,350,5]*0.04;
% settings.camtarget = [0,3,0];
settings.camtarget = [0,0,0];
settings.WriteOBJs = 1;
settings.MakeVideo = 1;
settings.SceneName = 'explorer1-a2';
settings.OBJDir = './objs/explorer1/';
settings.FramesToRecord = 7000;

td = simulate3D({mesha}, h, contactFinder, integrator, rigidificator, settings , energyModel);
save("explorer1AD_"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
writeTDcsv(td,"explorer1", ["_adaptive"]);
