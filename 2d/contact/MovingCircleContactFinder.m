classdef MovingCircleContactFinder < ContactFinder
    %CIRCLECONTACTFINDER Contact finder class that only computes contact
    %between meshes and a plane defined by normal and position
    
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
        function obj = MovingCircleContactFinder( radius, position, frictionCoefficient )
            obj.Radius = radius;
            obj.Center = position;
            if nargin >= 3
                obj.FrictionCoefficient = frictionCoefficient;
            end
            
%             obj.cfun = @(t) [ 0, 0.3*sin(t*6) ];
%             obj.dcdt = @(t) [ 0, 0.3*cos(t*6)*6 ];

%             obj.cfun = @(t) [ 0, 0.3*sin(t*6)^2 ];
%             obj.dcdt = @(t) [ 0, 0.3*2*sin(t*6)*cos(t*6)*6 ];

            obj.cfun = @(t) [ 0, 0.1*(cos(t*8)-1) * (mod( t*8, pi*2*10 ) < pi*2) ];
            obj.dcdt = @(t) [ 0, -0.1*sin(t*8)*8  * (mod( t*8, pi*2*10 ) < pi*2) ];

            obj.thetafun = @(t)  0.03*sin(t*6)    * (mod( t*6, pi*2*3 ) < pi*2);
            obj.dthetadt = @(t)  0.03*cos(t*6)*4  * (mod( t*6, pi*2*3 ) < pi*2);

            
        end
        
        function [J, phi, cInfo] = findContacts( obj, meshes, time )
            ps = vertcat(meshes.p);
            xs = ps(1:2:end);
            ys = ps(2:2:end);
            
            c = obj.Center + obj.cfun(time);
            r = obj.Radius;
            
            
            distance = sqrt((c(1) - xs).^2 + (c(2) - ys).^2);
            
            % indices of the ones under the plane
            idx = distance <= r;
            
            if ( all(idx == 0) )
                J = zeros(0,numel(ps));
                phi = [];
                cInfo = contactInfo.empty;
                return;
            end
            
            % constraint value (penetration)
            phi = distance(idx) - r;
            normal = [xs(idx)-c(1), ys(idx)-c(2)] ./ distance(idx);
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
                % will have other point based factors if drdt and dthetadt are not zero
                cInfo(i).velocity = -[dot(normal(i,:),obj.dcdt(time)); dot(tangent(i,:),obj.dcdt(time)) - obj.dthetadt(time)*r ]; 
            end         
        end
        
        function render( obj, time )
            c = obj.Center + obj.cfun(time);
            r = obj.Radius;
            theta = obj.thetafun(time);
            ct = cos(theta);
            st = sin(theta);
            R = [ct st; -st ct];
            p1 = r*R*[1;0] + c';
            p2 = r*R*[-1;0] + c';
            p3 = r*R*[0;1] + c';
            p4 = r*R*[0;-1] + c';
            th = linspace(0, pi*2,100);
            xdata = [c(1)+r*cos(th), nan, p1(1), p2(1), nan, p3(1), p4(1) ];
            ydata = [c(2)+r*sin(th), nan, p1(2), p2(2), nan, p3(2), p4(2) ];
            if ( obj.plotHandle == 0 )
                hold on;        
                obj.plotHandle = plot( xdata, ydata,  'Color','k','LineWidth', 0.5 );
            else
                obj.plotHandle.XData = xdata;
                obj.plotHandle.YData = ydata;
            end
        end
    end
end

