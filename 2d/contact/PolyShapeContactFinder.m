classdef PolyShapeContactFinder < ContactFinder
    % MESHCONTACTFINDER Class for finding contact between meshes.
    %   The shape can only have one boundary in the current implementation.
    properties
        finderPoly             % polyshape of this obstacle
        edgePoints             % points forming the polyshape
        minmaxM2               % AABB of this polyshape obstacle
        m2BoundaryEdgePoint1   
        m2BoundaryEdgePoint2
    end
    
    methods
        function obj = PolyShapeContactFinder( polyShapeNodes, frictionCoefficient )
            % PolyShapeContactFinder Construct an obstacle from the polygon
            % formed by the provided points.  Contacts with this obstacle
            % will be assigned the given frictionCoefficient.
            obj.finderPoly = polyshape( polyShapeNodes );
            obj.edgePoints = polyShapeNodes;
            obj.m2BoundaryEdgePoint1 = obj.edgePoints(1:1:end,:);
            obj.m2BoundaryEdgePoint2 = [obj.edgePoints(2:1:end,:);obj.edgePoints(1,:)];
            obj.minmaxM2 = minmax(obj.m2BoundaryEdgePoint1');
            if nargin >= 2
                obj.FrictionCoefficient = frictionCoefficient;
            end
        end
        
       function [Jc, phi, cInfo] = findContacts( obj, meshes, ~ )
            % Input k meshes N1,...,Nk DOFS
            % Output Jc will have N1+N2+...+Nk columns with contact * 2  rows
            % phi will be a column vector of the size contacts
            % contactInfo is a cell array with size equal to number of contacts
            
            mdofs = cell( 1, numel(meshes) );
            dofs = 0;
            for m = 1:numel(meshes)
                N = meshes(m).N*2;
                mdofs{m} = dofs+1:dofs+N;
                dofs = dofs + N;
            end
            
            phi = zeros( 0, 1 );
            Jc = sparse( 0, dofs );
            cInfo = contactInfo.empty;
            % find collisions
            for m1 = 1:numel(meshes)
                [Jc1,phiOut, cInfoOut] = findMeshContacts( obj, meshes(m1), mdofs{m1});                    
                if ( isempty(phiOut) )
                    continue;
                end
                Jcn = zeros( size(Jc1,1), dofs );
                Jcn( :, mdofs{m1} ) = Jc1;
                Jc = [Jc;Jcn];
                phi = [phi;phiOut];
                cInfo = [ cInfo, cInfoOut ];
            end
        end
        
        function [ Jc1, phi, cInfo ] = findMeshContacts( obj, mesh1, mdofsm1 )
            % Compute contact Jacobians for points on the boundary of mesh1
            % that fall within the boundary of mesh2.
            m1edges = vertcat( mesh1.boundaryEdges{:} );
            m1BoundaryPoints = [ mesh1.p(m1edges(:,1)*2-1), mesh1.p(m1edges(:,1)*2) ];
            
            % Collinear points will be removed 
            % great for efficiency, but produces a warning
            MSGID = 'MATLAB:polyshape:repairedBySimplify';
            warning('off', MSGID);
            
            m1InAABBm2 = ...
                m1BoundaryPoints(:,1) >= obj.minmaxM2(1,1) & ...
                m1BoundaryPoints(:,1) <= obj.minmaxM2(1,2) & ...
                m1BoundaryPoints(:,2) >= obj.minmaxM2(2,1) & ...
                m1BoundaryPoints(:,2) <= obj.minmaxM2(2,2);
            in = m1InAABBm2;
            if sum(in) <= 0
                Jc1 = [];
                Jc2 = [];
                phi = [];
                cInfo = contactInfo.empty;
                return;
            end
            in(m1InAABBm2) = isinterior( obj.finderPoly, m1BoundaryPoints(m1InAABBm2,1), m1BoundaryPoints(m1InAABBm2,2) );
            
            % if no interpenetration, exit early.
            if( all( in == 0 ) )
                Jc1 = [];
                phi = [];
                cInfo = {};
                return;
            end
            
            % Find phi as minimum distance, and choose opposite part of 
            % contact as that closest point on the edge of the other 
            % geometry.  This will violate equal and opposite friction 
            % forces, but is but probably OK if there isn't a lot of
            % interpenetration. 
            numContacts = sum(in~=0);
            m1ContactPoints = m1BoundaryPoints( in, : ); % contact points on bondary of mesh 1
            phi = inf( numContacts, 1 );
            m2ContactEdgeIndices = -ones( numContacts, 1 );  % index of corresponding edge on mesh 2

            % check all "in" contact points of shape1 with all boundary
            % edges of shape2.  This would be much better if it were 
            % vectorized overedges rather than points as we'll typically
            % have few "in" points and many edges. 
            for i = 1:size( obj.edgePoints, 1 )
                v1 = obj.m2BoundaryEdgePoint1( i, : );
                v2 = obj.m2BoundaryEdgePoint2( i, : );
                d = point_to_line_segment_distance( m1ContactPoints, v1, v2 );                
                % need to keep the minimum distance
                ix = find( d < phi );
                phi(ix) = d(ix);
                m2ContactEdgeIndices(ix) = i; % keep track of closest m1 edge index
            end

            % these interpenetraiton depths must be negative for
            % baumgarte to work!
            phi = -phi;
            
            % now for each in point, find the alpha, normal, tangent
            v1 = obj.m2BoundaryEdgePoint1( m2ContactEdgeIndices, : );
            v2 = obj.m2BoundaryEdgePoint2( m2ContactEdgeIndices, : );
            alphas = sum((m1ContactPoints-v1).*(v2-v1),2) ./ sum((v2-v1).*(v2-v1),2);
            alphas(alphas>1) = 1;
            alphas(alphas<0) = 0;
            m2ContactPoints = v1.*(1-alphas) + v2.*(alphas);  % contact point ons on mesh 2
            
            normals = m2ContactPoints - m1ContactPoints;
            normals = normals ./ sqrt(sum( normals.*normals, 2)); 
            tangents = [ - normals(:,2), normals(:,1) ];
            
            
            %init rows (velocity,friction)
            Jc1 = zeros( numContacts*2, 2*mesh1.N );
            Jc2 = zeros( numContacts*2, 2*size(obj.edgePoints,1) );

            
            % To build Jacobians, note that m1ContactPoints come directly 
            % from DOF indices via their index, while m1ContactPoints come
            % from the nodes of boundary edges.

            % Here, indices are of edge corresponding to mesh1 boundary edges
            % and recall the point was taken as the first index of the edge!
            m1EdgeIndices = find(in);
                        
            row = 1;
            cInfo = contactInfo.empty;
            for i = 1:numContacts
                %velocity
                Jc1( row, m1edges(m1EdgeIndices(i),1)*2 )   = normals(i,2);
                Jc1( row, m1edges(m1EdgeIndices(i),1)*2-1 ) = normals(i,1);      
                row = row + 1;
                
                %friction
                Jc1( row, m1edges(m1EdgeIndices(i),1)*2 )   = tangents(i,2);
                Jc1( row, m1edges(m1EdgeIndices(i),1)*2-1 ) = tangents(i,1);
                row = row + 1;
                
                % could sort these... but these will likewise always
                % be found in the same order
                indices = mdofsm1( m1edges(m1EdgeIndices(i),1) );
                cInfo(i) = contactInfo( ...
                    (m1ContactPoints(i,:) + m2ContactPoints(i,:))/2, ...
                    normals(i,:), ...
                    tangents(i,:), ...
                    obj.FrictionCoefficient, ...
                    indices, ...
                    obj.ID);
            end
            Jc1 = sparse( Jc1 );
        end
        
        function render( obj, ~ )
            if ( obj.plotHandle ~= 0 )
                return
            end
            hold on;
            obj.plotHandle = plot(obj.finderPoly,'FaceColor','green');
        end
    end
end