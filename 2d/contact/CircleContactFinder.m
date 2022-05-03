classdef CircleContactFinder < ContactFinder
    %CIRCLECONTACTFINDER Contact finder class that only computes contact
    %between meshes and a plane defined by normal and position
    
    properties
        Radius
        Center
    end
    
    methods
        function obj = CircleContactFinder( radius, position, frictionCoefficient )
            obj.Radius = radius;
            obj.Center = position;
            if nargin >= 3
                obj.FrictionCoefficient = frictionCoefficient;
            end
        end
        
        function [J, phi, cInfo] = findContacts( obj, meshes, time )
            ps = vertcat(meshes.p);
            xs = ps(1:2:end);
            ys = ps(2:2:end);
            
            distance = sqrt((obj.Center(1) - xs).^2 + (obj.Center(2) - ys).^2);
            
            % indices of the ones under the plane
            idx = distance <= obj.Radius;
            
            if ( all(idx == 0) )
                J = zeros(0,numel(ps));
                phi = [];
                cInfo = contactInfo.empty;
                return;
            end
            
            % constraint value (penetration)
            phi = distance(idx) - obj.Radius;
            normal = [xs(idx)-obj.Center(1), ys(idx)-obj.Center(2)] ./ distance(idx);
            tangent = [-normal(:,2), normal(:,1)];
            
            indices = find(idx);

            %add friction
            rown  = (1:2:2*numel(phi))';
            rowt  = (2:2:2*numel(phi))';
            colx = (indices*2-1);
            coly = (indices*2);
            J = sparse( ...
                [ rown;         rown;        rowt;         rowt  ], ...
                [ colx;         coly;        colx;         coly  ], ...
                [ normal(:,1);  normal(:,2); tangent(:,1); tangent(:,2) ], 2*numel(phi), numel(ps) );

                        
            cInfo = contactInfo.empty; %cell( numel(phi), 1 );
            for i = 1:numel(phi) 
                cInfo(i) = contactInfo( [xs(indices(i)), ys(indices(i))], normal(i,:), tangent(i,:), obj.FrictionCoefficient, indices(i), obj.ID);
            end         
        end
        
        function render( obj, ~ )
            if ( obj.plotHandle ~= 0 )
                return
            end
            theta = linspace(0,pi*2,100);
            c = obj.Center;
            r = obj.Radius;
            hold on;
            obj.plotHandle = plot( c(1)+r*cos(theta), c(2)+r*sin(theta),  'Color','k','LineWidth', 0.5 );
        end
    end
end

