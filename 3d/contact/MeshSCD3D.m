classdef MeshSCD3D < ContactFinder
    % MeshSCD3D Class for mesh self collision detection.
    %   The SCD test is not general, and instead works only for the case
    %   where the mesh consists of multiple connected components, and tests
    %   the boundaries betweeen those components.
    
    properties
        quality = false; % speed instead of quality by default
        skipset = [];
        epsilon = 0.0001;
    end
    
    methods
        function obj = MeshSCD3D( frictionCoefficient, quality , skipset)
            if nargin >= 1
                obj.FrictionCoefficient = frictionCoefficient;
            end
            if nargin >= 2
                obj.quality = quality;
            end
            if nargin >=3
                obj.skipset = skipset;
            end
        end
        
        function [Jc, phi, cInfo] = findContacts( obj, mesh,~)
            % Input mesh with N DOFS
            % Output Jc will have N columns with contact * 2  rows
            % phi will be a column vector of the size contacts
            % contactInfo is a cell array with size equal to number of contacts
            
            phi = zeros( 0, 1 );
            Jc = zeros( 0, mesh.N*3 );
            cInfo = contactInfo3D.empty;
            % find collisions
            for bes1 = 1:numel(mesh.facesSets)
                if any(bes1 == obj.skipset)
                    continue
                end
                others = bes1+1:numel(mesh.facesSets);
                if ( obj.quality )
                    others = 1:numel(mesh.facesSets);
                end
                for bes2 = others
                    if bes1 == bes2 || any(bes2 == obj.skipset)
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
            %Just to make sure that Jc is sparse and not mess up the 
            Jc = sparse(Jc); 
        end
        
        function render(obj,~) 
        end
        
        function [ Jc, phi, cInfo ] = findMeshContacts( obj, mesh, bes1, bes2 )
            % Compute contact Jacobians for points on the boundary of mesh1
            % that fall within the boundary of mesh2.
            
            bes1edges = mesh.facesSets{bes1};
            bes1Vert = 1:mesh.mergeMeshVerticeCount(bes1);
            bes2edges = mesh.facesSets{bes2};
            bes2Vert = 1:mesh.mergeMeshVerticeCount(bes2);
            
            adjustedBes1 = bes1Vert+ mesh.facesSetsOffset(bes1);
            adjustedBes2 = bes2Vert+ mesh.facesSetsOffset(bes2);
            
            s1BoundaryPoints     = [ mesh.p(adjustedBes1*3-2), mesh.p(adjustedBes1*3-1), mesh.p(adjustedBes1*3) ];
            s2BoundaryPoints     = [ mesh.p(adjustedBes2*3-2), mesh.p(adjustedBes2*3-1), mesh.p(adjustedBes2*3) ];

            %doing this is faster than minmax() for some reason
            AABB2 = [min( s2BoundaryPoints); max(s2BoundaryPoints)]'; % AABB 
            s1InAABB2 = ...
                s1BoundaryPoints(:,1) >= AABB2(1,1) & ...
                s1BoundaryPoints(:,1) <= AABB2(1,2) & ...
                s1BoundaryPoints(:,2) >= AABB2(2,1) & ...
                s1BoundaryPoints(:,2) <= AABB2(2,2) & ...
                s1BoundaryPoints(:,3) >= AABB2(3,1) & ...
                s1BoundaryPoints(:,3) <= AABB2(3,2) ;
            in = s1InAABB2;
            if sum(in) <= 0
                Jc = [];
                phi = [];
                cInfo = contactInfo3D.empty;
                return;
            end
            
            adjustedBes1= adjustedBes1(in);
            % Find phi as minimum distance, and choose opposite part of 
            % contact as that closest point on the edge of the other 
            % geometry.  This will violate equal and opposite friction 
            % forces, but is but probably OK if there isn't a lot of
            % interpenetration. 
            m1ContactPoints = s1BoundaryPoints( in, : ); % contact points on bondary of mesh 1

            [S,I,C,N] = signed_distance(m1ContactPoints, s2BoundaryPoints, bes2edges);

            % these interpenetraiton depths must be negative for
            % baumgarte to work!
            in = S < 0 & ~isnan(N(:,1)) & ~isnan(N(:,2)) & ~isnan(N(:,3));
            
            m1ContactPoints = m1ContactPoints(in,:);
            
            if sum(in) <= 0
                Jc = [];
                phi = [];
                cInfo = contactInfo3D.empty;
                return;
            end
            
            phi = S(in);
            phi(phi > 0) = 0;
            
            normals = N(in,:);
            tangents = zeros(size(normals));
            tangents2 = zeros(size(normals));
            for i = 1: size(normals,1)
                tangent = null(normals(i,:));
                tangents2(i,:) = tangent(:,2)';
                tangents(i,:) = tangent(:,1)';
            end
            numContacts = sum(in);
            
            %init rows (velocity,friction)
            Jc = zeros( numContacts*3, 3*mesh.N );
          
            % To build Jacobians, note that m1ContactPoints come directly 
            % from DOF indices via their index, while m1ContactPoints come
            % from the nodes of boundary edges.

            % Here, indices are of edge corresponding to mesh1 boundary edges
            % and recall the point was taken as the first index of the edge!
            m1ContactVertIndex = adjustedBes1(in);
            C = C(in,:);
            I = I(in);

            row = 1;
            cInfo = contactInfo3D.empty;
            
            for i = 1:numContacts
                verticesOnfaceHit = bes2edges(I(i),:);
                verticesAdjustedIds = adjustedBes2(verticesOnfaceHit);
                points = [mesh.p(verticesAdjustedIds*3-2),mesh.p(verticesAdjustedIds*3-1),mesh.p(verticesAdjustedIds*3)];
                dist = zeros(size(points));
                for j = 1:3
                    dist(:,j) = points(j,:) - C(i,:);
                end
                magnitude = sqrt(sum(dist.^2,2));
                
                alphas = 1 - magnitude/sum(magnitude);

                %velocity
                Jc( row, m1ContactVertIndex(i)*3 )    = normals(i,3);
                Jc( row, m1ContactVertIndex(i)*3-1 )  = normals(i,2);
                Jc( row, m1ContactVertIndex(i)*3-2 )  = normals(i,1);
                Jc( row, verticesAdjustedIds(1)*3 )   = -normals(i,3)' *alphas(1);
                Jc( row, verticesAdjustedIds(1)*3-1 ) = -normals(i,2)' *alphas(1);
                Jc( row, verticesAdjustedIds(1)*3-2 ) = -normals(i,1)' *alphas(1);
                Jc( row, verticesAdjustedIds(2)*3 )   = -normals(i,3)' *alphas(2);
                Jc( row, verticesAdjustedIds(2)*3-1 ) = -normals(i,2)' *alphas(2);
                Jc( row, verticesAdjustedIds(2)*3-2 ) = -normals(i,1)' *alphas(2);
                Jc( row, verticesAdjustedIds(3)*3 )   = -normals(i,3)' *alphas(3);
                Jc( row, verticesAdjustedIds(3)*3-1 ) = -normals(i,2)' *alphas(3);
                Jc( row, verticesAdjustedIds(3)*3-2 ) = -normals(i,1)' *alphas(3);
                row = row + 1;
                
                %friction
                Jc( row, m1ContactVertIndex(i)*3 )       = tangents(i,3);
                Jc( row, m1ContactVertIndex(i)*3-1 )     = tangents(i,2);
                Jc( row, m1ContactVertIndex(i)*3-2 )     = tangents(i,1);
                Jc( row, verticesAdjustedIds(1)*3 )   = -tangents(i,3)' *alphas(1);
                Jc( row, verticesAdjustedIds(1)*3-1 ) = -tangents(i,2)' *alphas(1);  
                Jc( row, verticesAdjustedIds(1)*3-2 ) = -tangents(i,1)' *alphas(1);     
                Jc( row, verticesAdjustedIds(2)*3 )   = -tangents(i,3)' *alphas(2);
                Jc( row, verticesAdjustedIds(2)*3-1 ) = -tangents(i,2)' *alphas(2);  
                Jc( row, verticesAdjustedIds(2)*3-2 ) = -tangents(i,1)' *alphas(2);  
                Jc( row, verticesAdjustedIds(3)*3 )   = -tangents(i,3)' *alphas(3);
                Jc( row, verticesAdjustedIds(3)*3-1 ) = -tangents(i,2)' *alphas(3);  
                Jc( row, verticesAdjustedIds(3)*3-2 ) = -tangents(i,1)' *alphas(3);  
                row = row + 1;
                
                                %friction
                Jc( row, m1ContactVertIndex(i)*3 )       = tangents2(i,3);
                Jc( row, m1ContactVertIndex(i)*3-1 )     = tangents2(i,2);
                Jc( row, m1ContactVertIndex(i)*3-2 )     = tangents2(i,1);
                Jc( row, verticesAdjustedIds(1)*3 )   = -tangents2(i,3)' *alphas(1);
                Jc( row, verticesAdjustedIds(1)*3-1 ) = -tangents2(i,2)' *alphas(1);  
                Jc( row, verticesAdjustedIds(1)*3-2 ) = -tangents2(i,1)' *alphas(1);     
                Jc( row, verticesAdjustedIds(2)*3 )   = -tangents2(i,3)' *alphas(2);
                Jc( row, verticesAdjustedIds(2)*3-1 ) = -tangents2(i,2)' *alphas(2);  
                Jc( row, verticesAdjustedIds(2)*3-2 ) = -tangents2(i,1)' *alphas(2);  
                Jc( row, verticesAdjustedIds(3)*3 )   = -tangents2(i,3)' *alphas(3);
                Jc( row, verticesAdjustedIds(3)*3-1 ) = -tangents2(i,2)' *alphas(3);  
                Jc( row, verticesAdjustedIds(3)*3-2 ) = -tangents2(i,1)' *alphas(3);  
                row = row + 1;
                
                % could sort these... but these will likewise always
                % be found in the same order
                indices = [m1ContactVertIndex(i), verticesAdjustedIds];
                cInfo(i) = contactInfo3D([(m1ContactPoints(i,:) + C(i,:))/2], ...
                    normals(i,:), ...
                    tangents(i,:), ...
                    obj.FrictionCoefficient, ...
                    indices, ...
                    obj.ID, ...
                    tangents2(i,:));
            end         
        end
    end
end