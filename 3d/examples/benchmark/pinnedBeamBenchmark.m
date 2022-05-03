    clear;
    close all;
    
testOutputFile = "resolutionTest.csv";
fileID = fopen(testOutputFile,'w');
nbytes = fprintf(fileID,'Test, Adaptive, Time\n',1);
fclose(fileID);


runTest(testOutputFile, "beam28",true);
% runTest(testOutputFile, "beam28",false);

runTest(testOutputFile,"beam730",true);
% runTest(testOutputFile,"beam730",false);
% 
runTest(testOutputFile,"beam913",true);
% runTest(testOutputFile,"beam913",false);
% 
runTest(testOutputFile,"beam1264",true);
% runTest(testOutputFile,"beam1264",false);
% 
runTest(testOutputFile,"beam2248",true);
% runTest(testOutputFile,"beam2248",false);
% 
runTest(testOutputFile,"beam5232",true);
% runTest(testOutputFile,"beam5232",false);
% 
runTest(testOutputFile,"beam8325",true);
% runTest(testOutputFile,"beam8325",false);
% 
runTest(testOutputFile,"beam16106",true);
runTest(testOutputFile,"beam16106",false);
% 
runTest(testOutputFile,"beam44608",true);
runTest(testOutputFile,"beam44608",false);
% 
runTest(testOutputFile,"beam116573",true);
runTest(testOutputFile,"beam116573",false);

function runTest(testFile, beamSizeName, isAdaptive)
    h = 0.01; % time step

    % MESHES

    rho = 10;
    nu = 0.45;      % Poisson ratio: close to 0.5 for rubber
    E = 1e6;     % Young's modulus: 0.01e9 approximate for rubber
    [mu, lambda] = toLame( nu, E );
    alpha0 = 0.01;   % Rayleigh factor on M
    alpha1 = 0.01;  % Rayleigh factor on K
    tMaterial = TriangleMaterial(rho, mu, lambda, alpha0, alpha1);

    mesh = meshLoader(beamSizeName, [], tMaterial);
    mesh.setRigidTransform( [90,90,0], [2,0,0]); % put it higher in the world frame
    % pinning tris
    xPos = mesh.p(1:3:end);
    pinInd = find(xPos == max(xPos));
    mesh.pin(sort(pinInd));

    if isAdaptive
        mesh = AdaptiveMesh3D(mesh);
    end

    rigidificator = EDot3DMexRigidificator();

    rigidificator.RigidificationThreshold = 2e-4;
    rigidificator.ElastificationThreshold = 5e-4; 
    rigidificator.FrameCount = 3;

    integrator = BackwardEuler3D();
    integrator.Compliance = 0.0;

    energyModel = StVenantKirchoff3DEnergy();

    planeContactFinder = NullContactFinder(3);
    contactFinder = {planeContactFinder};

    % CONFIG
    settings = Simulation3DSettings();
%     settings.PlotEDotHist = 1;
    settings.PlotSkip = 3000;%makes it headless
    settings.campos=[-30,40,5]*.5;
    settings.camtarget=[1,1,0];
    settings.FramesToRecord = 3000;

    td = simulate3D( mesh, h, contactFinder, integrator, rigidificator, settings, energyModel);
    
    textBool = "false";
    if isAdaptive
        textBool = "true";
    end
    
    fileID = fopen(testFile,'a');
    nbytes = fprintf(fileID,'%s, %s, %5d \n',beamSizeName, textBool, td.SimulationLoopFull);
    fclose(fileID);
    
    clear;
    close all;
end
