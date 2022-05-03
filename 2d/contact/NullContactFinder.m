classdef NullContactFinder < ContactFinder
    %Null equivalent for a contactFinder. This does not collision detection
    properties
        dimension
    end
    
    methods
        function obj = NullContactFinder(dimension)
            if nargin < 1
                dimension = 2;
            end
            obj.dimension = dimension;
        end
        function [J, phi, cInfo] = findContacts(obj, meshes, ~ )
            J = zeros(0, sum([meshes.N]) * obj.dimension);
            phi = zeros(0, 1);
            cInfo = {};
        end
    end
end

