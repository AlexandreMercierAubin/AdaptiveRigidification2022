classdef PlaneContactFinderMultiFriction < ContactFinder
    %PLANECONTACTFINDER Contact finder class that only computes contact
    %between meshes and a plane defined by normal and position
    
    properties
        Normal
        D
    end
    
    methods
        function obj = PlaneContactFinderMultiFriction(normal, position, frictionCoefficient)
            % PlaneContactFinder Takes normal, position, and mu to create a
            % new plane contact finder.
            obj.Normal = normal ./ norm(normal, 2);
            obj.D = dot(normal, position);
            
            if nargin >= 3
                obj.FrictionCoefficient = frictionCoefficient;
            end
        end
        
        function [J, phi, cInfo] = findContacts( obj, meshes, ~ )
            ps = vertcat(meshes.p);
            xs = ps(1:2:end);
            ys = ps(2:2:end);
            
            distance = xs * obj.Normal(1) + ys * obj.Normal(2) - obj.D;
            
            % indices of the ones under the plane
            idx = distance < 0;
            if ( all(idx == 0) )
                J = zeros(0,numel(ps));
                phi = [];
                cInfo = contactInfo.empty;
                return;
            end
            
            % constraint value (penetration)
            phi = distance(idx);
            normal = repmat( obj.Normal, numel(phi), 1 );
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
                if ( mod( floor(xs(indices(i))*2), 2 ) == 0 )
                    mu = 0;
                else
                    mu = obj.FrictionCoefficient;
                end
                cInfo(i) = contactInfo( [xs(indices(i)), ys(indices(i))], normal(i,:), tangent(i,:),mu, indices(i), obj.ID);
            end     
            
%                         norms = [idx, idx] .* obj.Normal;
% 
%             % for each contact constraint, we will need a friction constraint
%             tans = [idx, idx] .* tangent;
%             vals = permute(reshape([norms, tans]', 2, 2, []), [2 1 3]);
%             cells = num2cell(vals, [1 2]);
% 
%             
%             % build J matrix by attributing a row to each constraint
%             J = blkdiag(sparse(cells{1}), cells{2:end});
%             
%             % remove all rows that are only zeros 
%             J = J(any(J, 2), :);
%             
%             
%             indices = find(idx);
%             cInfo = contactInfo.empty;
%             for i = 1:numel(phi) 
%                 cInfo(i) = contactInfo( [xs(indices(i)), ys(indices(i))], obj.Normal, tangent, obj.FrictionCoefficient, indices(i), obj.ID);
%             end            
        end
        
        function render( obj, ~ )
            if ( obj.plotHandle ~= 0 )
                return
            end
            % plane equation is  x * obj.Normal(1) + y * obj.Normal(2) - obj.D = 0
            % so y = - N(1)/N(2) * x + D / N(2)
            m = -obj.Normal(1) / obj.Normal(2);
            b = obj.D / obj.Normal(2);
            x = -10:0.5:10;
            x = x';
            y = m*x + b;
            xd = [ x(1:2:end-1) x(2:2:end) nan(numel(x(2:2:end)),1) ];
            yd = [ y(1:2:end-1) y(2:2:end) nan(numel(y(2:2:end)),1) ];
            xd = reshape(xd',[],1);
            yd = reshape(yd',[],1);
            obj.plotHandle = plot( xd, yd,'b', 'LineWidth', 3);
            xd = [ x(2:2:end-1) x(3:2:end) nan(numel(x(3:2:end)),1) ];
            yd = [ y(2:2:end-1) y(3:2:end) nan(numel(y(3:2:end)),1) ];
            xd = reshape(xd',[],1);
            yd = reshape(yd',[],1);
            plot( xd, yd,'r', 'LineWidth', 3);
            %axis equal
        end
    end
end

