classdef ContactFinder < handle
    % CONTACTFINDER Finds contacts of the meshes in a simulation.
    %   Abstract class that specify a contract on finding contacts for
    %   meshes in the simulation. 
    
    properties
        FrictionCoefficient 
        ID
        plotHandle = 0
    end
    
    methods
        function obj = ContactFinder()
            persistent contactFinderCounter;
            if isempty(contactFinderCounter)
                contactFinderCounter = 1;
            else
                contactFinderCounter = contactFinderCounter + 1;
            end
            obj.ID = contactFinderCounter;
            
           obj.FrictionCoefficient = 0;
        end
        
        function render( obj, time )
            % RENDER renders this ContactFinder for plotting if needed, only
            % relevant for contact finders that introduce geometry like
            % plane and such.
            
            % nothing by default
        end
        
        function [F,V] = getObjPositionFaces(obj, time)
             % Allows the output of an obj file with the rendering or
             % approximate shape of the contact object
            
            % nothing by default
            F = [];
            V = [];
        end
    end
    
    methods (Abstract = true)
        [J, phi, cInfo] = findContacts( obj, meshes, time )
    end
end

