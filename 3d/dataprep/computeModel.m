function [V,F,U,q] = computeModalBasis(filename, varargin)
    % Input:
    %   filename  holding volume or surface mesh (in m)
    %             (.stl, .obj, .msh)
    %             if a volume mesh is provided, the corresponding surface 
    %             mesh with the same name should be in the same folder.
    %
    % Output: 
    %   .msh      generated volume mesh if a surface mesh is provided (generated by tetgen)
    %   .stl      computed surface mesh if a volume mesh is provided 
    %
    % Options:    (currently disabled)
    %   'Scale'   currently not functioning preoperly, default is 1
    %   'Draw'    followed by true or false, default is true
    %
    % Requirement: 
    % 1- Install tetgen: http://wias-berlin.de/software/tetgen/ 
    
    if nargin<1
        error('Filename required. Should be a surface mesh ".stl" or ".obj".');
    end
    
    %%%% OPTIONIAL PARAMS %%%%
    %% Default values
    
    scale = 1;
    draw = true;
    
    %% Map of parameter names to variable names
%     params = containers.Map( ...
%     {'NbModes','Boxes','Indices','Density','Young','Scale','Draw'},{'k','boxes','indices','density','young','scale','draw'});
%     v = 1;
%     while v <= numel(varargin)
%         param_name = varargin{v};
%         if isKey(params, param_name)
%               assert(v+1<=numel(varargin));
%               v = v+1;
%               % Trick: use feval on anonymous function to use assignin to this workspace 
%               feval(@()assignin('caller',params(param_name),varargin{v}));
%         else
%              error('Unsupported parameter: %s',varargin{v});
%         end
%         v=v+1;
%     end
    
    %% LOAD or, GENENERATE and EXPORT VOLUME MESH 
    [filepath,name,ext] = fileparts(filename);
    if ext == ".msh"
        [V,T,F] = readMSH(filename);
        V = V * scale;
        if isfile( filepath + "/" + name + ".obj" )
            [SV,SF] = loadMesh( filepath + "/" + name + ".obj" );
        elseif isfile(filepath+"/"+name+".stl")
            [SV,SF] = readSTL( filepath+"/"+name+".stl","JoinCorners",1);
        else
            error("Surface mesh with name "+name+" not found in "+filepath);
            return;
        end  
        SV = SV * scale;
    else
        if ext == ".stl"
            [SV,SF] = readSTL(filename,"JoinCorners",1); % Option is to remove duplicate vertices
        else
            [SV,SF] = loadMesh(filename);
        end
        SV = SV * scale;
        filename = filepath+"/"+name+".msh";
        [V,T,F] = tetgen(SV, SF);
        writeMSH(filename,V,T,F);
    end

    figure(1)
    clf;
    subplot(1,2,1);
    h = patch('vertices',V,'faces',F,'facecolor',[.5,.5,.5],'edgecolor',[0,0,0.9]);
    alpha 0.5;
    axis equal;
    axis off;
    lighting phong;
    camlight infinite;
    
    % This is super slow because of order correct transparency? 
    % (or might just be slow if there are tons of tetrahedra) 
    subplot(1,2,2);
    h2 = tetramesh(T,V, 'edgecolor','none');
    alpha 0.5;
    axis equal;
    axis off;
    lighting phong;
    camlight infinite;
    
    
end
