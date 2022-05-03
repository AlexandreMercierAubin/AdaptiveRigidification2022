    clear;
    close all;
    
testOutputFile = "sitffnessTest.csv";
fileID = fopen(testOutputFile,'w');
nbytes = fprintf(fileID,'Test, Adaptive, Time\n',1);
fclose(fileID);
runTest(testOutputFile,"beam2248",true,1e7);
runTest(testOutputFile,"beam2248",true,1e6);
runTest(testOutputFile,"beam2248",true,1e5);
runTest(testOutputFile,"beam2248",true,1e4);
runTest(testOutputFile,"beam2248",true,1e3);

testOutputFile = "RifidificationTest.csv";
fileID = fopen(testOutputFile,'w');
nbytes = fprintf(fileID,'Test, Adaptive, Time\n',1);
fclose(fileID);
runTestRigidification(testOutputFile,"beam2248",true,1e7);
runTestRigidification(testOutputFile,"beam2248",true,1e6);
runTestRigidification(testOutputFile,"beam2248",true,1e5);
runTestRigidification(testOutputFile,"beam2248",true,1e4);
runTestRigidification(testOutputFile,"beam2248",true,1e3);
runTestRigidification(testOutputFile,"beam2248",true,1e2);
runTestRigidification(testOutputFile,"beam2248",true,1e1);

testOutputFile = "ElastificationTest.csv";
fileID = fopen(testOutputFile,'w');
nbytes = fprintf(fileID,'Test, Adaptive, Time\n',1);
fclose(fileID);
runTestElastification(testOutputFile,"beam2248",true,1e7);
runTestElastification(testOutputFile,"beam2248",true,1e6);
runTestElastification(testOutputFile,"beam2248",true,1e5);
runTestElastification(testOutputFile,"beam2248",true,1e4);
runTestElastification(testOutputFile,"beam2248",true,1e3);
runTestElastification(testOutputFile,"beam2248",true,1e2);
runTestElastification(testOutputFile,"beam2248",true,1e1);

testOutputFile = "nuTest.csv";
fileID = fopen(testOutputFile,'w');
nbytes = fprintf(fileID,'Test, Adaptive, Time\n',1);
fclose(fileID);
runTestNu(testOutputFile,"beam2248",true,0.45);
runTestNu(testOutputFile,"beam2248",true,0.48);
runTestNu(testOutputFile,"beam2248",true,0.42);
runTestNu(testOutputFile,"beam2248",true,0.4);
runTestNu(testOutputFile,"beam2248",true,0.38);
runTestNu(testOutputFile,"beam2248",true,0.36);
runTestNu(testOutputFile,"beam2248",true,0.33);

testOutputFile = "alpha0Test.csv";
fileID = fopen(testOutputFile,'w');
nbytes = fprintf(fileID,'Test, Adaptive, Time\n',1);
fclose(fileID);
runTestAlpha0(testOutputFile,"beam2248",true,0.01);
runTestAlpha0(testOutputFile,"beam2248",true,0.1);
runTestAlpha0(testOutputFile,"beam2248",true,0.001);
runTestAlpha0(testOutputFile,"beam2248",true,0.0001);

testOutputFile = "alpha1Test.csv";
fileID = fopen(testOutputFile,'w');
nbytes = fprintf(fileID,'Test, Adaptive, Time\n',1);
fclose(fileID);
runTestAlpha1(testOutputFile,"beam2248",true,0.01);
runTestAlpha1(testOutputFile,"beam2248",true,0.1);
runTestAlpha1(testOutputFile,"beam2248",true,0.001);
runTestAlpha1(testOutputFile,"beam2248",true,0.0001);

function runTest(testFile, beamSizeName, isAdaptive, E)
    h = 0.01; % time step

    % MESHES

    rho = 10;
    nu = 0.45;      % Poisson ratio: close to 0.5 for rubber
    
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
    nbytes = fprintf(fileID,'%d, %s, %5d \n',E, textBool, td.SimulationLoopFull);
    fclose(fileID);
    
    clear;
    close all;
end

function runTestRigidification(testFile, beamSizeName, isAdaptive, R)
    h = 0.01; % time step

    % MESHES

    rho = 10;
    nu = 0.45;      % Poisson ratio: close to 0.5 for rubber
    E = 1e6;
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

    rigidificator.RigidificationThreshold = R;
    rigidificator.ElastificationThreshold = 5e-4; 
    rigidificator.FrameCount = 3;

    integrator = BackwardEuler3D();
    integrator.Compliance = 0.0;

    energyModel = StVenantKirchoff3DEnergy();

    planeContactFinder = NullContactFinder(3);
    contactFinder = {planeContactFinder};

    % CONFIG
    settings = SimulationSettings();
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
    nbytes = fprintf(fileID,'%d, %s, %5d \n',R, textBool, td.SimulationLoopFull);
    fclose(fileID);
    
    clear;
    close all;
end

function runTestElastification(testFile, beamSizeName, isAdaptive, El)
    h = 0.01; % time step

    % MESHES

    rho = 10;
    nu = 0.45;      % Poisson ratio: close to 0.5 for rubber
    E = 1e6;
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
    rigidificator.ElastificationThreshold = El; 
    rigidificator.FrameCount = 3;

    integrator = BackwardEuler3D();
    integrator.Compliance = 0.0;

    energyModel = StVenantKirchoff3DEnergy();

    planeContactFinder = NullContactFinder(3);
    contactFinder = {planeContactFinder};

    % CONFIG
    settings = SimulationSettings();
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
    nbytes = fprintf(fileID,'%d, %s, %5d \n',El, textBool, td.SimulationLoopFull);
    fclose(fileID);
    
    clear;
    close all;
end

function runTestNu(testFile, beamSizeName, isAdaptive, nu)
    h = 0.01; % time step

    % MESHES
    rho = 10;
    E = 1e6;
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
    settings = SimulationSettings();
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
    nbytes = fprintf(fileID,'%d, %s, %5d \n',nu, textBool, td.SimulationLoopFull);
    fclose(fileID);
    
    clear;
    close all;
end

function runTestAlpha0(testFile, beamSizeName, isAdaptive, alpha0)
    h = 0.01; % time step

    % MESHES
    rho = 10;
    E = 1e6;
    nu = 0.45;
    [mu, lambda] = toLame( nu, E );
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
    settings = SimulationSettings();
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
    nbytes = fprintf(fileID,'%d, %s, %5d \n',alpha0, textBool, td.SimulationLoopFull);
    fclose(fileID);
    
    clear;
    close all;
end

function runTestAlpha1(testFile, beamSizeName, isAdaptive, alpha1)
    h = 0.01; % time step

    % MESHES
    rho = 10;
    E = 1e6;
    nu = 0.45;
    [mu, lambda] = toLame( nu, E );
    alpha0 = 0.01;
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
    settings = SimulationSettings();
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
    nbytes = fprintf(fileID,'%d, %s, %5d \n',alpha1, textBool, td.SimulationLoopFull);
    fclose(fileID);
    
    clear;
    close all;
end

