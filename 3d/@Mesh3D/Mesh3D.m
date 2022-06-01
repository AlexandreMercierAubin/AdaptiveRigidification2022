classdef Mesh3D < handle
    properties
        N   % number of DOFs
        p   % node positions N by 3
        pSurface %node positions of the mesh skin
        SurfaceWeights %weights of the tets on the mesh skin, set as [] to deactivate the feature
        Surface2Tet % array of size #surface vertices storing corresponding tet for coloring the vertices | must be either used for all meshes or none
        SurfaceFaces
        surfaceVertexBaseColor %color of individual surface vertices (pre-rigidification)
        prevp
        p0  % initial node positions
        v   % node velocities
        v0  % initial node velocities
        f   % force accumulator for nodes
        
        t      % tet (for making elements)
        tint32
        el     % elements
        faces  % cell array of boundaries (for collision detection)
        % Every mesh has at least one set of boundary edges
        % Each edge is a line, first column is the index of the first point
        % second column is the index of the second point
        % NOTE: order is important!!!  Treating edges as directed, the mesh
        % is to the left andt the outside is to the right, e.g., CCW order
        % for convex polygons.
        faceEdges % vertex ids of each edge as rows edges*2 [vertex1, vertex2]. No duplicate. 
        triEdges % edge IDS of a tri. tri by 3
        facesSets % which sets of boundaries can be checked with each other for collision
        facesSetsOffset
        facesTriList
        objectDOFs %cell array of lists that containts the dofs of each individual objects
        objectDOFsPinned %same, but with pinned DOFs removed
        mergeMeshVerticeCount
        kdRestEdgeAverageOppositeHeight % bend stiffness * edgeRestLength / restEdgeAverageOppositeHeight

        cacheMD5 % MD5 solely use to detect change in merged meshes scene creation
        cacheMergeMaterials %solely use to detect change in merged meshes scene creation
        
        dphidx

        Graph                   % Tet adjacency graph (across edges)
        AdjagencyMatrix         % Adjacency matrix used to build the graph (not used otherwise?)
        TetsPerParticle         % cellarray of the neighboring tets of each individual vertex
        
        pinned       % flags for pinned vertex indices
        pinnedInds   % indices of pinned vertices
        pinnedDOFs   % pinned indices in the 3 * N state vector
        unpinnedDOFs % unpinned indices in the 3 * N state vector
        pinnedTets   % pinned tets (tets with at least a pinned particle)
        isTetPinned  % flags for pinned tets 
        stablePinnedElements %tris with two pins or tets with 3 pins
        animationDOFs = [];
        animationInds = [];
        
        activeDOFs
        ActiveBRows % indices of rows in B which are elastic (active)
        
        
        mass    % diag of mass matrix (likely more convenient?)
        M       % sparse mass matrix, 3 * N by 3 * N
        Mii
        Minv
        alpha0mass     % diag of lumped Rayleigh alpha0 multiplied with lumped mass
        Md             % sparse matrix version of alpha0mass     
        Mdii
        
        B              % kinematic relationship between vertices and element deformation gradients F
        Bii
        Bt
        
        initialB
        lagrangeMults  % lambdas for the compliant feedback constraint version of elasticity

        %Cf   % symbolic d2psid2F function - deprectated use mCSTVK
        elV         % quick per element volume access for C computation
        elMu        % quick per element Lame access for C computation
        elLambda    % quick per element lambda access for C computation
        elAlpha1    % quick per element alpha 1 acces
        
        bigAlpha1   % sparse matrix for combining with full C to build Kd
        
        RenderPatch        % plot of the elements
        FaceLines      % plot of the boundaries
        
        materials           % material per tet
        materialIndex       % tet index to material
        faceMaterialIndex
        
        valence %valence of the vertices
        lineHandles %debugging line handles kept for a more efficient plotting of lines
    end
    
    methods
        function obj = Mesh3D(p, t, attributes, materials, boundaryFaces, boundaryTIndices, SurfaceWeights, SurfaceFaces, pSurface, Surface2Tet)
            % MESH Constructs a mesh from given points and tets.
            %   Prepares for simulation with provided materials.
            %
            %   Note that if you use meshes saved in a .mat file, you will
            %   need to regenerate the mesh if new fields get added
            %   
            if nargin == 2 || (nargin >= 3 && size(attributes,2) <= 0)
                attributes = ones(size(t,1),1);
            end
            
            if nargin <= 3 || (nargin >= 4 && size(materials,2) <= 0)
                defaultMaterial = TriangleMaterial();
            else
                defaultMaterial = materials;
            end
            
            if nargin >= 2 && size(attributes,2) > size(materials,2)
                ids = setdiff(1:size(attributes,2),1:size(materials,2));
                defaultMaterial(ids) = TriangleMaterial();
            end
            
            if isa(p, 'Mesh3D')
               % copy constructor
                fns = properties(p);
                for i = 1:numel(fns)
                    obj.(fns{i}) = p.(fns{i});
                end
                return;
            end
            
            if nargin < 5
                [boundaryFaces,obj.facesSets, boundaryTIndices] = makeBoundaries3D(p, t);
            end
            
            if nargin < 7
                SurfaceWeights = [];
                SurfaceFaces = [];
                pSurface = [];
                Surface2Tet = [];
            end
            
            N = size(p, 1);
            obj.N = N;                 % number of nodes
            obj.objectDOFs = {1:N*3};
            obj.objectDOFsPinned = obj.objectDOFs;
            obj.p = reshape(p', N * 3, 1); % position state
            obj.prevp = obj.p;
            obj.p0 = obj.p;           % initial positions
            obj.v = zeros(N * 3, 1);      % velocity state
            obj.v0 = zeros(N * 3, 1);     % initial velocities
            obj.f = zeros(N * 3, 1);      % force accumulator
            obj.mergeMeshVerticeCount = [obj.N];
            
            obj.SurfaceWeights = SurfaceWeights;
            obj.SurfaceFaces = SurfaceFaces;
            obj.pSurface = pSurface;
            obj.Surface2Tet = Surface2Tet;
            obj.t = t;                 % tetrahedrons
            obj.tint32 = int32(t);

            obj.el = obj.makeElements3D(p, t);
            if size(obj.t,2) == 3
                obj.referenceSpaceNormal = zeros(9*size(obj.t,1),size(obj.t,1)*3);
                for i = 1:numel(obj.el)
                    n = obj.el(i).referenceSpaceNormal;
                    obj.referenceSpaceNormal(i*9-8:i*9-6, i*3-2) = n;
                    obj.referenceSpaceNormal(i*9-5:i*9-3,i*3-1) = n;
                    obj.referenceSpaceNormal(i*9-2:i*9,i*3) = n;
                end
                obj.referenceSpaceNormal = sparse(obj.referenceSpaceNormal);
            end
            
            obj.facesTriList = boundaryTIndices;
            obj.faces = boundaryFaces;
            obj.facesSets = {boundaryFaces};
            obj.facesSetsOffset = [0];
            
            duplicateEdge = containers.Map('KeyType','char','ValueType','double');
            obj.faceEdges = [];
            for i = 1: size(obj.faces,1)
                vert1 = obj.faces(i,1);
                vert2 = obj.faces(i,2);
                vert3 = obj.faces(i,3);
                if ~duplicateEdge.isKey(char([vert2,vert3])) && ~duplicateEdge.isKey(char([vert3,vert2])) 
                    obj.faceEdges = [obj.faceEdges; vert2, vert3];
                    duplicateEdge(char([vert3,vert2])) = true;
                end
                
                if ~duplicateEdge.isKey(char([vert1,vert2])) && ~duplicateEdge.isKey(char([vert2,vert1])) 
                    obj.faceEdges = [obj.faceEdges; vert1, vert2];
                    duplicateEdge(char([vert1,vert2])) = true;
                end
                
                if ~duplicateEdge.isKey(char([vert1,vert3])) && ~duplicateEdge.isKey(char([vert3,vert1])) 
                    obj.faceEdges = [obj.faceEdges; vert1, vert3];
                    duplicateEdge(char([vert1,vert3])) = true;
                end
            end
            
            obj.triEdges = zeros(size(obj.t,1),3);
            if size(obj.t,2) == 3
                for i = 1:size(obj.faceEdges,1)
                    edges = obj.faceEdges(i,:);
                    for j = 1:size(obj.t,1)
                        edgeInTri = any(edges(1)==obj.t(j,:)) && any(edges(2)==obj.t(j,:));
                        if edgeInTri
                            for edgeSlot = 1:3
                                if obj.triEdges(j,edgeSlot) == 0
                                    obj.triEdges(j,edgeSlot) = i;
                                    break;
                                end
                            end
                        end
                    end
                end
            end
            
            %computing valence
            obj.valence = zeros(obj.N,1);
            for i = 1:size(obj.t,1)
                for j = 1:size(obj.t,2)
                    vertex = obj.t(i,j);
                    obj.valence(vertex) = obj.valence(vertex) + 1;
                end
            end
            
            obj.AdjagencyMatrix = sparse(size(t, 1));
            obj.TetsPerParticle = cell(N, 1, 1);
            
            for i = 1:1:size(t, 1)
                build(i);
            end
            
            % could not find a better way to vectorize this but if there
            % is, its probably faster.
            function build(i)
                if size(obj.t,2) == 3
                    inds = ismember(1:N, t(i, :));
                    obj.TetsPerParticle(inds) = cellfun(@(x) [x, i], obj.TetsPerParticle(inds), 'UniformOutput', false);
                    % Check all triangles with higher index, to see which
                    % indices are shared with our triangle... and if exactly
                    % two are shared then we share an edge with that triangle!
                    % That is, we've found the adjacent triangles.
                    matching = find(sum(ismember(t(i + 1:size(t, 1), :), t(i, :)),2) == 2) + i;
                
                else
                    inds = ismember(1:N, t(i, :));

                    obj.TetsPerParticle(inds) = cellfun(@(x) [x, i], obj.TetsPerParticle(inds), 'UniformOutput', false);

                    matching = find(sum(ismember(t(i + 1:size(t, 1), :), t(i, :)),2) == 3) + i;
                end
                obj.AdjagencyMatrix(i, matching) = 1;
                obj.AdjagencyMatrix(matching, i) = 1; 
            end
            
            obj.Graph = graph(obj.AdjagencyMatrix);
            
            obj.pinned = zeros(N, 1);   % flags pinned indices
            obj.pinnedInds = [];        % list of pinned node IDs
            obj.pinnedDOFs = [];
            obj.unpinnedDOFs  = 1:obj.N*3;
            obj.pinnedTets = [];
            obj.isTetPinned = zeros(size(obj.t,1),1);
            obj.activeDOFs = 1:3 * N;
            obj.ActiveBRows = 1:size(t, 1) * 9;

            V = reshape(obj.p, 3, size(obj.p, 1) / 3)';
            
            if size(obj.t,2) == 3
                obj.B = computeB3D(V, obj.t);
                obj.Bii = obj.B(:,obj.unpinnedDOFs);
                obj.Bt = sparse(obj.B');
                obj.initialB = obj.B;
                obj.dphidx = linear_tri3dmesh_dphi_dX(V, obj.t);
            else
                obj.B = computeB3D(V, obj.t);
                obj.Bt = obj.B';
                obj.Bii = obj.B(:,obj.unpinnedDOFs);
            end
            obj.lagrangeMults = zeros(size(obj.B, 1), 1);
            obj.elV = volume(V,obj.t);

            obj.updateMaterials(attributes,[defaultMaterial]);
            
            alpha1s = reshape( repmat( obj.elAlpha1, 9, 1 ), [], 1 );
            obj.bigAlpha1 = sparse( 1:numel(alpha1s), 1:numel(alpha1s), alpha1s );
            
            obj.RenderPatch = 0;
            obj.FaceLines = {};
        end
        
        function mergeMesh(obj, mesh) 
            % MergeMesh adds a given mesh to this mesh
            % this is not necessarily efficient in that it will recompute
            % things that were computed before, but it all precomputation
            % so perhaps we don't need to worry so much.
            
            N1 = obj.N;
            N2 = mesh.N;
            numTets1 = size(obj.t,1);
            
            obj.N = N1 + N2;  
            obj.objectDOFs = {obj.objectDOFs{:}, mesh.objectDOFs{:} + (N1*3)};
            obj.objectDOFsPinned = obj.objectDOFs;
            obj.p = [ obj.p; mesh.p ];
            obj.prevp = obj.p;
            obj.p0 = [ obj.p0; mesh.p0 ];
            obj.v = [ obj.v; mesh.v ];
            obj.v0 = [ obj.v0; mesh.v0 ];
            obj.f = [ obj.f; mesh.f ];
            obj.mergeMeshVerticeCount = [obj.mergeMeshVerticeCount, mesh.mergeMeshVerticeCount];
            % triangle indices of the added mesh are offset by the number
            % of particles in this mesh
            obj.t = [ obj.t; mesh.t + N1 ];  % tets
            obj.Surface2Tet = [obj.Surface2Tet, mesh.Surface2Tet + numTets1];
            obj.tint32 = int32(obj.t);
            
            obj.SurfaceWeights = sparse(blkdiag(obj.SurfaceWeights,mesh.SurfaceWeights));
            nFaces = numel(obj.pSurface)/3;
            obj.SurfaceFaces = [obj.SurfaceFaces;nFaces+mesh.SurfaceFaces];
            obj.pSurface = [obj.pSurface;mesh.pSurface];
            
            obj.valence = [obj.valence; mesh.valence];
            obj.el = obj.makeElements3D( reshape( obj.p0', 3, [] )', obj.t );
   
            numBoundarySets1 = numel(obj.facesSets);
            for j = 1:numel(mesh.facesSets)
                obj.facesSets{numBoundarySets1+j} = mesh.facesSets{j};
                obj.facesSetsOffset(numBoundarySets1+j) = obj.mergeMeshVerticeCount(numBoundarySets1+j-1) + obj.facesSetsOffset(numBoundarySets1+j-1);
            end
            obj.faces = [obj.faces;mesh.faces + N1]; 
            obj.facesTriList = [obj.facesTriList;mesh.facesTriList + numTets1]; 

            obj.AdjagencyMatrix = blkdiag(obj.AdjagencyMatrix,mesh.AdjagencyMatrix);
            obj.TetsPerParticle = [obj.TetsPerParticle;mesh.TetsPerParticle];
                       
            obj.activeDOFs = 1:3 * obj.N;
            obj.ActiveBRows = 1:size(obj.t, 1) * 9;
            
            % silly to recompute, but whatever... 
            V = reshape(obj.p, 3, size(obj.p, 1) / 3)';
            obj.B = blkdiag(obj.B,mesh.B);
            obj.Bt = obj.B';
            obj.Bii = obj.B(:,obj.unpinnedDOFs);

            obj.lagrangeMults = zeros(size(obj.B, 1), 1);
            obj.elV = volume(V,obj.t);
  
            numMat1 = numel(obj.materials);
            obj.updateMaterials([ obj.materialIndex; mesh.materialIndex + numMat1 ], [ obj.materials, mesh.materials ])
            
            alpha1s = reshape( repmat( obj.elAlpha1, 9, 1 ), [], 1 );
            obj.bigAlpha1 = sparse( 1:numel(alpha1s), 1:numel(alpha1s), alpha1s );

            obj.RenderPatch = 0;
            obj.FaceLines = {};
            obj.pin(mesh.pinnedInds + N1);
        end
        
        function clone = clone(obj)
            clone = Mesh3D(obj);
        end
        
        function updateMaterials( obj, newAttributes, newMaterials )
            % UPDATEMATERIALS Allows for new materials to be assigned to a
            % mesh after loading.
            obj.materials = newMaterials;
            obj.materialIndex = newAttributes;
            obj.faceMaterialIndex = obj.materialIndex(obj.facesTriList); 
            obj.elMu = [ obj.materials(obj.materialIndex(:)).mu ];
            obj.elLambda = [ obj.materials(obj.materialIndex(:)).lambda ];
            obj.elAlpha1 = [ obj.materials(obj.materialIndex(:)).alpha1 ];
            updateMass(obj);
            obj.M = sparse( 1:3*obj.N, 1:3*obj.N, obj.mass ); % sparse matrix form for convenience
            obj.Mii = obj.M(obj.unpinnedDOFs,obj.unpinnedDOFs);
            obj.Minv = sparse( 1:3*obj.N, 1:3*obj.N, obj.mass.^-1 );
            updateAlpha0(obj);
            obj.Md = sparse( 1:3*obj.N, 1:3*obj.N, obj.alpha0mass); % sparse matrix form for convenience
            obj.Mdii = obj.Md(obj.unpinnedDOFs,obj.unpinnedDOFs);
            
            obj.surfaceVertexBaseColor = [];
            if ~isempty(obj.SurfaceWeights)
                obj.surfaceVertexBaseColor = zeros(size(obj.Surface2Tet,1),3); 
                for i=1:numel(obj.materialIndex(obj.Surface2Tet))
                    tetID = obj.Surface2Tet(i);
                    matID = obj.materialIndex(tetID);
                    obj.surfaceVertexBaseColor(i,:) = obj.materials(matID).color;
                end
            end
        end
            
        function [B, gamma] = getB(obj, ~)
           B = obj.B;
           gamma = speye(obj.N * 3);
        end
        
        function M = getM(mesh)
            M = mesh.M;
        end
        
        function f = getCurrentForce(mesh)
            f = mesh.f;
        end
        
        function v = getCurrentVelocity(mesh)
            v = mesh.v;
        end
        
        function linearMomentum = getLinearMomentum( mesh )
            linearMomentum = sum( reshape(mesh.mass.*mesh.v, 3, mesh.N), 2 );
        end
        
        function labelTriangles( mesh )
             pt = reshape( mesh.p, 3, mesh.N );
             for i=1:numel(mesh.el)
                 pos = sum( pt( :, mesh.el(i).t ), 3 ) / 4;
                 text( pos(1), pos(2), pos(3), ""+i );
             end        
        end

        function labelVertices( mesh )
             pt = reshape( mesh.p, 3, mesh.N );
             for i=1:mesh.N
                 pos = pt(:,i);
                 text( pos(1), pos(2), pos(3), ""+i );
             end        
        end

        function updateMass(mesh)
            % UPDATEMASS is used to contsruct the mass damping and is
            % typically only called on construction of the mesh
            mesh.mass = computeMass(mesh);
        end
        
        function updateAlpha0( mesh )
            % UPDATEALPHA0 is used to contsruct the mass damping and is
            % typically only called on construction of the mesh
            
            alpha0tet = [mesh.materials(mesh.materialIndex).alpha0];
            %values are declared per tet, but we want them per node
            alpha0vector = zeros( mesh.N*3, 1 );
            for i = 1:size(mesh.t,1)
                for node = 1:size(mesh.t,2)
                    id = mesh.t(i,node);
                    value = alpha0vector(id*3) + 1/4 * alpha0tet(i);
                    alpha0vector(id*3) = value;
                    alpha0vector(id*3-1) = value;
                    alpha0vector(id*3-2) = value;
                end
            end
            mesh.alpha0mass = alpha0vector .* mesh.mass ;
        end
        
        %% public function prototypes
        prepare( obj, infMassForPinned )
        pin( obj, pinnedInds )
        updateParticles( obj, h, deltav )
        applyAcceleration( obj, acc )
        resetForce( obj )
        setRigidTransform( obj, degrees, translation, scale )
        setRigidMotion( obj, omega, radPerSec, velocity  )
        setRigidRotationFromMatrix(obj, rotationMatrix, forceCenterOfMass)
        computeActiveDOFs(obj) 
    end
    methods(Static)
        [W, T, SV, V, V2Htet] = voxelizer(V,FACES,side);
        [W, V2Htet] = meshSkinning(Vsurface,TetMesh3D);
        [el] = makeElements3D(V, T);
    end
end

