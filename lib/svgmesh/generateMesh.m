function mesh = generateMesh( filename, varargin )
    %GENERATEMESH Creates a triangle mesh from the given file, or drawing.
    %
    % From a given poly or a given svg:
    %     triangulates from contour
    % No input:
    %     same but will start with a drawing tool (svg file will be saved too)
    %
    % Input:
    %   filename  .svg or .poly file, or name of output without file extension
    % Output:
    %   .1.POLY .1.NODE and .1.ELE file that can be loaded with loadMeshFromPOLY 
    %
    % Options:
    %   'Draw'    followed by true or false to draw the mesh, default false
    %   'MaxArea' followed by maximum area constraint, default 0.01
    %
    % Examples:
    %   generateMesh( "data/model")                     % will launch a drawing tool and save the contour in model.svg
    %   generateMesh( "data/model.svg" )                % WARNING: an automatic scale of mesh assumes the mesh is in mm
    %   generateMesh( "data/model.svg", "Draw", true )  % triangulate and draw the mesh 
    %   generateMesh( "data/model.poly" )               % triangulates the given POLY file (e.g., hand made files)
    %
    % Dependencies: 
    %   https://github.com/alecjacobson/gptoolbox
    %   Image Processing Toolbox, if called without file extension (i.e., for drawing)
    %

    % Map of parameter names to variable names
    draw = false;
    maxArea = 0.01;
    minimumAngle = 30;
    params_to_variables = containers.Map( {'Draw', 'MaxArea', 'MinAngle'}, {'draw', 'maxArea', 'minAngle'} );
    v = 1;
    while v <= numel(varargin)
        param_name = varargin{v};
        if isKey(params_to_variables,param_name)
            assert(v+1<=numel(varargin));
            v = v+1;
            feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
        else
            error('Unsupported parameter: %s',varargin{v});
        end
        v=v+1;
    end
        
    [filepath,modelname,ext] = fileparts(filename);
    if modelname == ""
        modelname="model";
    end
    if filepath == ""
        filepath = ".";
    end
    
    points = [];
    edges = [];
    holes = [];
    triangulate = true;
    exportSVG = false;
    useProvidedPoly = false;

    if (ext == ".obj" || ext == ".stl")
        [points, edges] = load_mesh(filename);
        points(:,3) = []; edges(:,3) = []; 
    elseif (ext == ".svg")
        displayPlot = true;
        R = loadSVG(filename, 0.5, displayPlot);
        points = unique(R{1},'rows','stable');
        points = points * 1e-4;
    elseif (ext == ".poly")
        useProvidedPoly = true;
    elseif (ext == "")
        figure(1);
        axis([-2 2 -2 2]);
        axis equal;
        P = drawpolygon;
        points = P.Position;
        exportSVG = true;
    else
        error('Unknown mesh format: %s',ext);
    end
        
    if isempty(edges)
        N = size(points,1);
        s = 1:N;
        segments = [s;mod(s,N)+1];
        edges = segments';
    end
    
    if exportSVG
        saveSVG(filepath + modelname + ".svg", points, edges, [1,1,1], [0,0,0], 1);
    end
    
    if triangulate
        if (useProvidedPoly == true ) 
            system( pwd + "/lib/svgmesh/triangle.exe -p -q30 -A -a" + maxArea  + " " + filename );
            [TV,I] = readNODE(filepath + "/" + modelname + ".1.node");
            [TF,A] = readELE(filepath + "/" + modelname + ".1.ele");
        else 
            writePOLY_triangle( filepath + "/" + modelname + ".poly", points, edges, holes);
            system( pwd + "/lib/svgmesh/triangle.exe  -p -q" + minimumAngle + " -A -a" + maxArea  + " " + filepath + "/" + modelname + ".poly" );
            [TV,I] = readNODE( filepath + "/" + modelname + ".1.node" );
            [TF,A] = readELE( filepath + "/" + modelname + ".1.ele");
            %delete generateMeshTemp.*
        end
        
        if isempty(TF)
            TF = [];
        else
            % Triangle likes to use 1-indexed though .ele reader is 0-indexed
            if(( min(TF(:)) > 1) && (max(TF(:)) > size(TV,1)))
                TF = TF-1;
            elseif min(TF(:)) == 0 && max(TF(:)) < size(TV,1)
                TF = TF+1;
            end
        end
    end
    
    %% Example on how to use this to make a mesh...
    % actually... just use the loadMeshFromPOLY
    
%     % Make mesh with mechanical parameters
%     mu = young / 2 / (1 - nu);                  % shear modulus
%     lambda = young * nu / (1+nu) / (1-2*nu);    % Lambda parameter
%     alpha = 0.001;                          % Rayleigh mass
%     beta = 0.001;                           % Rayleigh stiffness
%     
%     mesh = Mesh( TV, TF, density, mu, lambda, alpha, beta );
%     save("data/" + modelname + ".MAT", "mesh");
     
    if draw
        clf;
        if numel(A) > 0
            colours = [ 0.5 0.5 0.5; 1 1 0.5; 0.5 1 1; 1 0.5 1; 1 0.5 0.5; 0.5 0.5 1; 0.5 1 0.5 ];
            faceColor = colours(mod(A,7)+1,:);
            patch('vertices', TV, 'faces', TF, 'edgecol', 'k',  'facecol', 'flat', 'FaceVertexCData', faceColor, 'FaceAlpha', .5, 'EdgeAlpha', .9 );
        else
            patch('vertices', TV, 'faces', TF, 'edgecol', 'k',  'facecol', [0.5 0.5 0.5], 'FaceAlpha', .5, 'EdgeAlpha', .9 );
        
        end
        axis equal
    end

end 



