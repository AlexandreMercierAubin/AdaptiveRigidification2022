classdef PlaneContactFinder3D < ContactFinder
    %PLANECONTACTFINDER Contact finder class that only computes contact
    %between meshes and a plane defined by normal and position
    properties
        Normal
        D
        position
        Tangent
        Tangent2
        planeSize = 1000;
		renderType = RenderType.wireframe;
    end
    
    methods
        function obj = PlaneContactFinder3D(normal, position, frictionCoefficient)
            % PlaneContactFinder Takes normal, position, and mu to create a
            % new plane contact finder.
            obj@ContactFinder();
            obj.Normal = normal / norm(normal);
            obj.D = dot(normal, position);
            obj.position = position;
            obj.Tangent = null(obj.Normal);
            obj.Tangent2 = obj.Tangent(:,2)';
            obj.Tangent = obj.Tangent(:,1)';
            obj.Tangent = obj.Tangent/norm(obj.Tangent);
            obj.Tangent2 = obj.Tangent2/norm(obj.Tangent2);
            
            if nargin >= 3
                obj.FrictionCoefficient = frictionCoefficient;
            end
        end
        
        function [J, phi, cInfo] = findContacts( obj, meshes, ~ )
            ps = vertcat(meshes.p);
            xs = ps(1:3:end);
            ys = ps(2:3:end);
            zs = ps(3:3:end);
            
            distance = (xs * obj.Normal(1) + ys * obj.Normal(2) + zs * obj.Normal(3)) - obj.D;
            
            % indices of the ones under the plane
            idx = distance <= 0;
            if ( all(idx == 0) )
                J = zeros(0,numel(ps));
                phi = [];
                cInfo = contactInfo3D.empty;
                return;
            end
            
            % constraint value (penetration)
            phi = distance(idx);
            normal = repmat( obj.Normal, numel(phi), 1 );
            tangent = repmat( obj.Tangent, numel(phi), 1 );
            tangent2 = repmat( obj.Tangent2, numel(phi), 1 );
            
            indices = find(idx);
            
            %add friction
            rown = (1:3:3*numel(phi))';
            rowt = (2:3:3*numel(phi))';
            rowt2 = (3:3:3*numel(phi))';
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
			
			if obj.renderType == RenderType.square
				[V,F] = obj.getObjPositionFaces(0);
				obj.plotHandle = patch('Faces',F,'Vertices',V,'FaceColor','red','EdgeColor', 'None');
				return;
			end
            nlines = 40;
            s = 10;
            tl1 = repmat(obj.position,nlines*3,1);
            tl1(1:3:end,:) = tl1(1:3:end,:) + repmat(obj.Tangent,nlines,1) .* linspace(-s,s,nlines)' + repmat(obj.Tangent2,nlines,1)*-s;
            tl1(2:3:end,:) = tl1(2:3:end,:) + repmat(obj.Tangent,nlines,1) .* linspace(-s,s,nlines)' + repmat(obj.Tangent2,nlines,1)*s;
            tl1(3:3:end,:) = nan; % breaks the connection in the lines
            tl2 = repmat(obj.position,nlines*3,1);
            tl2(1:3:end,:) = tl2(1:3:end,:) + repmat(obj.Tangent2,nlines,1) .* linspace(-s,s,nlines)' + repmat(obj.Tangent,nlines,1)*-s;
            tl2(2:3:end,:) = tl2(2:3:end,:) + repmat(obj.Tangent2,nlines,1) .* linspace(-s,s,nlines)' + repmat(obj.Tangent,nlines,1)*s;
            tl2(3:3:end,:) = nan; % breaks the connection in the lines
            tl = [ tl1; tl2 ];
            obj.plotHandle = plot3( tl(:,1), tl(:,2), tl(:,3), 'r' );

        end
        
        function [V,F] = getObjPositionFaces(obj, time)
            V(1,1:3) = obj.position - obj.Tangent*obj.planeSize - obj.Tangent2*obj.planeSize;
            V(2,1:3) = obj.position - obj.Tangent*obj.planeSize + obj.Tangent2*obj.planeSize;
            V(3,1:3) = obj.position + obj.Tangent*obj.planeSize - obj.Tangent2*obj.planeSize;
            V(4,1:3) = obj.position + obj.Tangent*obj.planeSize + obj.Tangent2*obj.planeSize;
            F(1,1:3) = [1,3,4];
            F(2,1:3) = [1,4,2];
        end
    end
end

