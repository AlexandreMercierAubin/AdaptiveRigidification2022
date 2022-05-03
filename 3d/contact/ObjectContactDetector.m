classdef ObjectContactDetector < ContactFinder
    %ObjectContactDetector collision detector that allows the collision
    %detection between meshes and an imported obj mesh
    
    properties
        V %vertices
        F %faces
        AABB
    end
    
    methods
        function obj = ObjectContactDetector(filename, degrees, scale, position, frictionCoefficient )
            [V,F,UV,TF,N,NF] = readOBJ(filename);
            R = RotationMatrixDegrees(degrees);
            
            obj.V = V;
            obj.V = (R*obj.V')';
            obj.V = (diag(scale)*obj.V')';
            obj.V(:,1) = obj.V(:,1) + position(1);
            obj.V(:,2) = obj.V(:,2) + position(2);
            obj.V(:,3) = obj.V(:,3) + position(3);
            
            obj.AABB = minmax( obj.V' );
            obj.F = F;
            obj.FrictionCoefficient = frictionCoefficient;
        end
        
        function [J, phi, cInfo] = findContacts( obj, meshes, time )
            ps = vertcat(meshes.p);
            xs = ps(1:3:end);
            ys = ps(2:3:end);
            zs = ps(3:3:end);
            
            % AABB
            s1InAABB2 = ...
                xs >= obj.AABB(1,1) & ...
                xs <= obj.AABB(1,2) & ...
                ys >= obj.AABB(2,1) & ...
                ys <= obj.AABB(2,2) & ...
                zs >= obj.AABB(3,1) & ...
                zs <= obj.AABB(3,2) ;
            idx = s1InAABB2;

            if ( sum(idx) == 0 )
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

