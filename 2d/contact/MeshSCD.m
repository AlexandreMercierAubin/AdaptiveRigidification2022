classdef MeshSCD < ContactFinder
    % MESHSCD Class for mesh self collision detection.
    %   The SCD test is not general, and instead works only for the case
    %   where the mesh consists of multiple connected components, and tests
    %   the boundaries betweeen those components.
    
    properties
        quality = false; % speed instead of quality by default
        epsilon = 1e-8;
    end
    
    methods
        function obj = MeshSCD( frictionCoefficient, quality )
            if nargin >= 1
                obj.FrictionCoefficient = frictionCoefficient;
            end
            if nargin >= 2
                obj.quality = quality;
            end
        end
        
        function [Jc, phi, cInfo] = findContacts( obj, mesh, ~)
            % Input mesh with N DOFS
            % Output Jc will have N columns with contact * 2  rows
            % phi will be a column vector of the size contacts
            % contactInfo is a cell array with size equal to number of contacts
            
            phi = zeros( 0, 1 );
            Jc = zeros( 0, mesh.N*2 );
            cInfo = contactInfo.empty;
            % find collisions
            for bes1 = 1:numel(mesh.boundaryEdgeSets)
                others = bes1+1:numel(mesh.boundaryEdgeSets);
                if ( obj.quality )
                    others = 1:numel(mesh.boundaryEdgeSets);
                end
                for bes2 = others  
                    if bes1 == bes2 %|| bes2 < bes1%tmp second condition
                        continue;
                    end
                    [JcOut, phiOut, cInfoOut] = findMeshContacts( obj, mesh, bes1, bes2 );
                    %meshes(m1), meshes(m2), mdofs{m1}, mdofs{m2} );                    
                    if ( isempty(phiOut) )
                        continue;
                    end
                    Jc = [ Jc; JcOut ];
                    phi = [ phi; phiOut ];
                    cInfo = [ cInfo, cInfoOut ];
                end
            end
            %Just to make sure that Jc is sparse and not mess up the solver later 
            Jc = sparse(Jc); 
        end
        
        function render(obj) 
        end
        
        function [ Jc, phi, cInfo ] = findMeshContacts( obj, mesh, bes1, bes2 )
            % Compute contact Jacobians for points on the boundary of mesh1
            % that fall within the boundary of mesh2.
            
            bes1edges = vertcat( mesh.boundaryEdges{ mesh.boundaryEdgeSets{bes1} } );
            bes2edges = vertcat( mesh.boundaryEdges{ mesh.boundaryEdgeSets{bes2} } );
            
            s1BoundaryPoints     = [ mesh.p(bes1edges(:,1)*2-1), mesh.p(bes1edges(:,1)*2) ];
            s2BoundaryPoints     = [ mesh.p(bes2edges(:,1)*2-1), mesh.p(bes2edges(:,1)*2) ];

            AABB2 = minmax( s2BoundaryPoints' ); % AABB
            s1InAABB2 = ...
                s1BoundaryPoints(:,1) >= AABB2(1,1) & ...
                s1BoundaryPoints(:,1) <= AABB2(1,2) & ...
                s1BoundaryPoints(:,2) >= AABB2(2,1) & ...
                s1BoundaryPoints(:,2) <= AABB2(2,2);
            in = s1InAABB2;
            if sum(in) <= 0
                Jc = [];
                phi = [];
                cInfo = contactInfo.empty;
                return;
            end
            
            % Collinear points will be removed 
            % great for efficiency, but produces a warning
            MSGID = 'MATLAB:polyshape:repairedBySimplify';
            warning('off', MSGID);
            mbes2 = mesh.boundaryEdgeSets{ bes2 };
            b2first = mesh.boundaryEdges{ mbes2(1) };
            shape2 = polyshape( mesh.p(b2first(:,1)*2-1), mesh.p(b2first(:,1)*2) );            
            for j=2:numel(mbes2)  
                edges = mesh.boundaryEdges{ mbes2(j) };
                assert( numel(mbes2) == 2 );
                shape2j = polyshape( mesh.p(edges(:,1)*2-1), mesh.p(edges(:,1)*2) );
                shape2 = subtract(shape2, shape2j);
            end        

            in(s1InAABB2) = isinterior( shape2, s1BoundaryPoints(s1InAABB2,1), s1BoundaryPoints(s1InAABB2,2) );
            
            % if no interpenetration, another opportunity to exit early.
            if( all( in == 0 ) )
                Jc = [];
                phi = [];
                cInfo = contactInfo.empty;
                return;
            end
            
            % Find phi as minimum distance, and choose opposite part of 
            % contact as that closest point on the edge of the other 
            % geometry.  This will violate equal and opposite friction 
            % forces, but is but probably OK if there isn't a lot of
            % interpenetration. 
            numContacts = sum(in~=0);
            m1ContactPoints = s1BoundaryPoints( in, : ); % contact points on bondary of mesh 1
            phi = inf( numContacts, 1 );
            m2ContactVertIndex1 = -ones( numContacts, 1 );  % index of first vertex of edge on mesh 2
            m2ContactVertIndex2 = -ones( numContacts, 1 );  % index of second vertex of edge on mesh 2
            

            % check all "in" contact points of boundary edge set 1 with all
            % boundary edges of boundary edge set 2.  This would be much
            % better if it were vectorized overedges rather than points as
            % we'll typically have few "in" points and many edges.
            for j = 1:numel(mbes2)
                edges = mesh.boundaryEdges{ mbes2(j) };
                for i = 1:size( edges, 1 )
                    v1 = [ mesh.p(edges(i,1)*2-1), mesh.p(edges(i,1)*2) ];
                    v2 = [ mesh.p(edges(i,2)*2-1), mesh.p(edges(i,2)*2) ];                    
                    d = point_to_line_segment_distance( m1ContactPoints, v1, v2 );                
                    % need to keep the minimum distance
                    ix = find( d < phi );
                    phi(ix) = d(ix);
                    m2ContactVertIndex1(ix) = edges(i,1);
                    m2ContactVertIndex2(ix) = edges(i,2);
                end
            end
            %polyshape does not filter exact 0 dist contacts so we do
            isZeroInterpenetration = phi < obj.epsilon;
            phi(isZeroInterpenetration) = [];
            m2ContactVertIndex1(isZeroInterpenetration) = [];
            m2ContactVertIndex2(isZeroInterpenetration) = [];

            % these interpenetraiton depths must be negative for
            % baumgarte to work
            phi = -phi;
            
            % now for each in point, find the alpha, normal, tangent
            v1 = [ mesh.p(m2ContactVertIndex1*2-1), mesh.p(m2ContactVertIndex1*2) ];
            v2 = [ mesh.p(m2ContactVertIndex2*2-1), mesh.p(m2ContactVertIndex2*2) ];
            alphas = (sum((m1ContactPoints-v1).*(v2-v1),2) ./ sum((v2-v1).*(v2-v1),2));
%             assert(~any(alphas>1));
%             assert(~any(alphas<0));
%             if any(alphas>1)
%                 test= 1;
%             end
%             if any(alphas<0)
%                 test = 2;
%             end
            m2ContactPoints = v1.*(1-alphas) + v2.*(alphas);  % contact point ons on mesh 2
            
            contactVectors = m2ContactPoints - m1ContactPoints;
            norms = vecnorm(contactVectors,2,2);

            normals = contactVectors ./ norms; 
            tangents = [ - normals(:,2), normals(:,1) ];

            %init rows (velocity,friction)
            Jc = zeros( numContacts*2, 2*mesh.N );
          
            % To build Jacobians, note that m1ContactPoints come directly 
            % from DOF indices via their index, while m1ContactPoints come
            % from the nodes of boundary edges.

            % Here, indices are of edge corresponding to mesh1 boundary edges
            % and recall the point was taken as the first index of the edge!
            m1ContactVertIndex = bes1edges(in,1);
            indices = [ m1ContactVertIndex, m2ContactVertIndex1, m2ContactVertIndex2 ];
            
            row = 1;
            cInfo = contactInfo.empty;
            for i = 1:numContacts
                %velocity
                Jc( row, m1ContactVertIndex(i)*2 )    = normals(i,2);
                Jc( row, m1ContactVertIndex(i)*2-1 )  = normals(i,1);
                Jc( row, m2ContactVertIndex1(i)*2 )   = -normals(i,2)'*(alphas(i));
                Jc( row, m2ContactVertIndex1(i)*2-1 ) = -normals(i,1)'*(alphas(i));
                Jc( row, m2ContactVertIndex2(i)*2 )   = -normals(i,2)'*(1-alphas(i));
                Jc( row, m2ContactVertIndex2(i)*2-1 ) = -normals(i,1)'*(1-alphas(i));                
                row = row + 1;
                
                %friction
                Jc( row, m1ContactVertIndex(i)*2 )       = tangents(i,2);
                Jc( row, m1ContactVertIndex(i)*2-1 )     = tangents(i,1);
                Jc( row, m2ContactVertIndex1(i)*2 )   = -tangents(i,2)'*(alphas(i));
                Jc( row, m2ContactVertIndex1(i)*2-1 ) = -tangents(i,1)'*(alphas(i));     
                Jc( row, m2ContactVertIndex2(i)*2 )   = -tangents(i,2)'*(1-alphas(i));
                Jc( row, m2ContactVertIndex2(i)*2-1 ) = -tangents(i,1)'*(1-alphas(i));
                row = row + 1;
                
                % could sort these... but these will likewise always
                % be found in the same order
                
                cInfo(i) = contactInfo( ...
                    mesh.p(m1ContactVertIndex(i)*2-1:m1ContactVertIndex(i)*2)', ...
                    normals(i,:), ...
                    tangents(i,:), ...
                    obj.FrictionCoefficient, ...
                    indices(i,:), ...
                    obj.ID);
                cInfo(i).pointAlpha = mesh.p(m2ContactVertIndex1(i)*2-1:m2ContactVertIndex1(i)*2);
                cInfo(i).pointAlphaInv = mesh.p(m2ContactVertIndex2(i)*2-1:m2ContactVertIndex2(i)*2);
                cInfo(i).normalAlpha = -normals(i,:).*(alphas(i));
                cInfo(i).normalAlphaInv = -normals(i,:).*(1-alphas(i));
                cInfo(i).tangentAlpha = -tangents(i,:).*(alphas(i));
                cInfo(i).tangentAlphaInv = -tangents(i,:).*(1-alphas(i));
            end         
        end
    end
end