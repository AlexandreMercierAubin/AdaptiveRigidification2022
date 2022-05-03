classdef FrameData
    % Data for a single frame used to present some debug information visually
    
    properties
        Contacts
        ContactNodePositions
        EDots
        FrameDuration
        Lambdas
        Velocities
        DeltaV
        p
        elasticBodiesCount
        rigidBodiesCount
        kinetic
        elasticKinetic
        rigidKinetic
    end
    
    methods
        function obj = FrameData()
            obj.Contacts = [];
            obj.ContactNodePositions = [];
            obj.EDots = [];
            obj.FrameDuration = 0.01;
            obj.Lambdas = [];
            obj.Velocities = [];
            obj.p = [];
            obj.kinetic = [];
            obj.elasticKinetic = [];
            obj.rigidKinetic = [];
        end
    end
end

