classdef SphereContactFinder < ContactFinder
    %SphereContactFinder Contact finder class that only computes contact
    %between meshes and a sphere defined by normal and position
    
    properties
        Radius
        Center
        plotRatio
        V
        F
    end
    
    methods
        function obj = SphereContactFinder( radius, position, frictionCoefficient )
            obj.Radius = radius;
            obj.Center = position;
            if nargin >= 3
                obj.FrictionCoefficient = frictionCoefficient;
            end
            obj.plotRatio = 1.0;
            [obj.V,obj.F] = subdivided_sphere(3);
        end
        
        function [J, phi, cInfo] = findContacts( obj, meshes, time )
            ps = vertcat(meshes.p);
            xs = ps(1:3:end);
            ys = ps(2:3:end);
            zs = ps(3:3:end);
            
            distance = sqrt((obj.Center(1) - xs).^2 + (obj.Center(2) - ys).^2 + (obj.Center(3) - zs).^2);
            
            % indices of the ones under the plane
            idx = distance <= obj.Radius;
            
            if ( all(idx == 0) )
                J = zeros(0,numel(ps));
                phi = [];
                cInfo = contactInfo3D.empty;
                return;
            end
            
            % constraint value (penetration)
            phi = (distance(idx) - obj.Radius);
            normal = [xs(idx) - obj.Center(1), ys(idx) - obj.Center(2), zs(idx) - obj.Center(3)] ./ distance(idx);
            tangent = zeros(size(normal));
            tangent2 = zeros(size(normal));
            for i = 1: size(normal,1)
                tangents = null(normal(i,:));
                tangent(i,:) = tangents(:,1)';
                tangent2(i,:) = tangents(:,2)';
            end
            
            
            indices = find(idx);

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
            V = obj.Center + obj.V*obj.Radius*obj.plotRatio;
            obj.plotHandle = patch('Faces',obj.F,'Vertices',V,'FaceColor','Blue');
        end
        
        function [V,F] = getObjPositionFaces(obj, time)
            F = obj.F;
            V = obj.Center + obj.V*obj.Radius*obj.plotRatio;
        end
    end
end

