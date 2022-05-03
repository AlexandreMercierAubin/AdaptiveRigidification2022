classdef PlatformContactDetector < ContactFinder
    %PlatformContactDetector collision detector that allows the collision
    %detection between meshes and a rectangular prism platform
    
    properties
        HalfWidth
        HalfHeight
        HalfDepth
        DirectionX
        DirectionY
        DirectionZ
        Center
        V %vertices
        F %faces
    end
    
    methods
        function obj = PlatformContactDetector( degrees, scale, position, frictionCoefficient )
            obj.HalfWidth = scale(1)/2;
            obj.HalfHeight = scale(2)/2;
            obj.HalfDepth = scale(3)/2;
            
            obj.Center = position;
            
            R = RotationMatrixDegrees(degrees);
            
            obj.DirectionX = (R*[1,0,0]')';
            obj.DirectionY = (R*[0,1,0]')';
            obj.DirectionZ = (R*[0,0,1]')';
            if nargin >= 8
                obj.FrictionCoefficient = frictionCoefficient;
            end
            
            obj.V(1,1:3) = obj.Center - obj.DirectionX * obj.HalfWidth - obj.DirectionY * obj.HalfHeight - obj.DirectionZ * obj.HalfDepth;
            obj.V(2,1:3) = obj.Center - obj.DirectionX * obj.HalfWidth - obj.DirectionY * obj.HalfHeight + obj.DirectionZ * obj.HalfDepth;
            obj.V(3,1:3) = obj.Center + obj.DirectionX * obj.HalfWidth - obj.DirectionY * obj.HalfHeight - obj.DirectionZ * obj.HalfDepth;
            obj.V(4,1:3) = obj.Center + obj.DirectionX * obj.HalfWidth - obj.DirectionY * obj.HalfHeight + obj.DirectionZ * obj.HalfDepth;
            obj.V(5,1:3) = obj.Center - obj.DirectionX * obj.HalfWidth + obj.DirectionY * obj.HalfHeight - obj.DirectionZ * obj.HalfDepth;
            obj.V(6,1:3) = obj.Center - obj.DirectionX * obj.HalfWidth + obj.DirectionY * obj.HalfHeight + obj.DirectionZ * obj.HalfDepth;
            obj.V(7,1:3) = obj.Center + obj.DirectionX * obj.HalfWidth + obj.DirectionY * obj.HalfHeight - obj.DirectionZ * obj.HalfDepth;
            obj.V(8,1:3) = obj.Center + obj.DirectionX * obj.HalfWidth + obj.DirectionY * obj.HalfHeight + obj.DirectionZ * obj.HalfDepth;
            
            obj.F = [1,6,5;
                     1,2,6;
                     2,8,6;
                     2,4,8;
                     4,8,7;
                     4,7,3;
                     3,7,5;
                     3,5,1;
                     6,8,7;
                     6,7,5;
                     2,4,3;
                     2,3,1]; 
        end
        
        function [J, phi, cInfo] = findContacts( obj, meshes, time )
            ps = vertcat(meshes.p);
            xs = ps(1:3:end);
            ys = ps(2:3:end);
            zs = ps(3:3:end);
            
            direction = [xs-obj.Center(1),ys-obj.Center(2),zs-obj.Center(3)];
            
            % indices of the points within the box
            idx = false(size(direction,1),1);
            for i = 1:size(direction,1)
                idx(i) = abs(dot(direction(i,:), obj.DirectionX)) <= obj.HalfWidth && ...
                        abs(dot(direction(i,:), obj.DirectionY)) <= obj.HalfHeight && ...
                        abs(dot(direction(i,:), obj.DirectionZ)) <= obj.HalfDepth;
            end
            if ( all(idx == 0) )
                J = zeros(0,numel(ps));
                phi = [];
                cInfo = contactInfo3D.empty;
                return;
            end
            
            points = [xs(idx),ys(idx),zs(idx)];
            
           [S,I,C,N] = signed_distance(points, obj.V, obj.F);

            % these interpenetraiton depths must be negative for
            % baumgarte to work!
            in = S < 0 & ~isnan(N(:,1)) & ~isnan(N(:,2)) & ~isnan(N(:,3));
            
            points = points(in,:);
            
            if sum(in) <= 0
                J = zeros(0,numel(ps));
                phi = [];
                cInfo = contactInfo3D.empty;
                return;
            end
            
            phi = S(in);
            phi(phi > 0) = 0;
            
            normal = N(in,:);
            tangent = zeros(size(normal));
            tangent2 = zeros(size(normal));
            for i = 1: size(normal,1)
                tangents = null(normal(i,:));
                tangent(i,:) = tangents(:,1)';
                tangent2(i,:) = tangents(:,2)';
            end
            
            indices = find(idx);
            indices = indices(in);

            %add friction
            rown  = (1:3:3*numel(phi))';
            rowt  = (2:3:3*numel(phi))';
            rowt2  = (3:3:3*numel(phi))';
            colx = (indices*3-2);
            coly = (indices*3-1);
            colz = (indices*3);
            J = sparse( ...
                [ rown;         rown;      rown;       rowt;         rowt;     rowt;   rowt2;         rowt2;     rowt2], ...
                [ colx;         coly;      colz;       colx;         coly;     colz;   colx;         coly;     colz;], ...
                [ normal(:,1);  normal(:,2); normal(:,3); tangent(:,1); tangent(:,2) ; tangent(:,3) ; tangent2(:,1); tangent2(:,2) ; tangent2(:,3)], 3*numel(phi), numel(ps) );

            cInfo = contactInfo3D.empty; %cell( numel(phi), 1 );
            for i = 1:numel(phi) 
                cInfo(i) = contactInfo3D( [xs(indices(i)), ys(indices(i)), zs(indices(i))], normal(i,:), tangent(i,:), obj.FrictionCoefficient, indices(i), obj.ID, tangent2(i,:));
            end             
        end
        
        function render( obj, ~ )
            if ( obj.plotHandle ~= 0 )
                return
            end
            
            hold on;
            obj.plotHandle = patch('Faces',obj.F,'Vertices',obj.V,'FaceColor','red');
        end
        
        function [V,F] = getObjPositionFaces(obj, time)
            F = obj.F;
            V = obj.V;
        end
    end
end

