classdef MovingSphereContactFinder < ContactFinder
    %WIP SphereContactFinder Contact finder class that only computes contact
    %between meshes and a sphere defined by normal and position. The sphere
    %moves according to a callback function
    
    properties
        Radius
        Center
        
        % NOTE: trajectory could be a property rather than hard coded...
        % could be a function that could be evaluated.
        cfun
        dcdt
        thetafun
        dthetadt
    end
    
    methods
        function obj = MovingSphereContactFinder( radius, position, frictionCoefficient )
            obj.Radius = radius;
            obj.Center = position;
            if nargin >= 3
                obj.FrictionCoefficient = frictionCoefficient;
            end
            
            obj.cfun = @(t) [ 0, 0.1*(cos(t*8)-1) * (mod( t*8, pi*2*10 ) < pi*2), 0];
            obj.dcdt = @(t) [ 0, -0.1*sin(t*8)*8  * (mod( t*8, pi*2*10 ) < pi*2), 0];

            obj.thetafun = @(t)  0.03*sin(t*6)    * (mod( t*6, pi*2*3 ) < pi*2);
            obj.dthetadt = @(t)  0.03*cos(t*6)*4  * (mod( t*6, pi*2*3 ) < pi*2);
        end
        
        function [J, phi, cInfo] = findContacts( obj, meshes, time )
            ps = vertcat(meshes.p);
            xs = ps(1:3:end);
            ys = ps(2:3:end);
            zs = ps(3:3:end);
            
            c = obj.Center + obj.cfun(time);
            distance = sqrt((c - xs).^2 + (c - ys).^2 + (c - zs).^2);
            
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
            normal = [xs(idx) - c(1), ys(idx) - c(2), zs(idx) - c(3)] ./ distance(idx);
            tangent = [-normal(:,2), normal(:,1), normal(:,3)];
            tangent2 = cross(tangent(:,:),normal(:,:),2);
            
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
                cInfo(i).velocity = -[dot(normal(i,:),obj.dcdt(time)); dot(tangent(i,:),obj.dcdt(time)) - obj.dthetadt(time)*distance; dot(tangent2(i,:),obj.dcdt(time)) - obj.dthetadt(time)*distance]; 
            end             
        end
        
        function render( obj, time )
            c = obj.Center + obj.cfun(time);
            r = obj.Radius*0.97;
            
            [x,y,z] = sphere;
            x = x * r;
            y = y * r;
            z = z * r;
            
            if ( obj.plotHandle == 0 )
                hold on;
                obj.plotHandle = surf(c(1)+x,c(2)+y,c(3)+z,'EdgeColor','none');
            else
                obj.plotHandle.XData = c(1)+x;
                obj.plotHandle.YData = c(2)+y;
                obj.plotHandle.ZData = c(3)+z;
            end
            
        end
    end
end

