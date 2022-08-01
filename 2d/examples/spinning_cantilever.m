clear;

h = 0.01; %time step

rho = 40;
nu = 0.4; % Poisson ratio: close to 0.5 for rubber
E = 2e4; % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.1; % Rayleigh factor on M
alpha1 = 0.1; % Rayleigh factor on K
material = TriangleMaterial( rho, mu, lambda, alpha0, alpha1 );

settings = SimulationSettings();
% settings.DrawTimings = 1;
settings.SceneName = "spinningCantilever";

resetMesh = true;
scale = [1,1];
rot = 0;
mesh1 = AdaptiveMesh(fetchPoly2D('cantileverP2',resetMesh, material, scale, rot, settings)); %twoTriSym, barP2, SixTri

mesh1.setRigidMotion( 3, [0,0] ); 
mesh2 = AdaptiveMesh(mesh1);

be = LDLBackwardEuler();
be.Gravity = 0;

rigid = EDiffMexRigidificator();
rigid.RigidificationThreshold = 1e-4;
rigid.ElastificationThreshold = 1e-3;
rigid.Preconditionner = "ICHOL";

% run a comparison
td = simulate( {mesh1,mesh2}, be, h, settings, rigid );

if ( numel(td) > 1 ) 
    figure(4); clf;
    td{1}.plotCompare( td{2}, 'Rng', 2 );  
    title('Simulation Time');
    figure(5); clf;
    td{1}.plotCompare( td{2} );
    title('Simulation Breakdown');
else 
    figure(4); clf;
    td{1}.plot( 'Rng', 2  );  
    title('Simulation Time');
    figure(5); clf;
    td{1}.plot();
    title('Simulation Breakdown');
end
