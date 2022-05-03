classdef Integrator < handle
    properties
        Gravity
        InfMassForPinned
        Name
        CustomForce 
        Baumgarte
        Compliance
        prevCInfo
        useFullAinv = false;
        PGSquicksolve = false;
        separateQuicksolveGravity = true; %might want this to be off in scenes with pinned dofs
    end
    
    methods
        function obj = Integrator()
            %INTEGRATOR Creates an integrator, object used to simulate
            %meshes
            obj.Gravity = -9.8;
            obj.InfMassForPinned = 0;
            obj.Name = 'Unamed';
            % specify vector or scalar if you want force added to your meshes at every frame
            obj.CustomForce = 0; 
            obj.Baumgarte =  10; %default baumgarte
            obj.Compliance = 0; %default compliance
            obj.prevCInfo = {};
            obj.useFullAinv = false;
        end
        
        function setComplianceAndBaumgarteFromKandB( obj, h, spring, damping )
             %computation of compliance according to the step size (1/h^2) http://www.ode.org/ode-latest-userguide.html#sec_3_8_2
            ERP = (h*spring)/(h*spring+damping);
            CFM = 1.0/(h*spring + damping);
            obj.setComplianceAndBaumgarteFromERPandCFM(h, ERP, CFM);
        end
        
        function setComplianceAndBaumgarteFromERPandCFM( obj,h, ERP, CFM)
            obj.Compliance = CFM;
            obj.Baumgarte = ERP/h;
        end
    end
    
    methods (Abstract = true)
        integrate( obj, meshes, h, Jc, phi, cInfo, settings, cache, td, animationScripts, frame );
        %INTEGRATE integrates the specified system
        %
        %   meshes: list of meshes
        %   h: time step
        %   Jc: contact constraints jacobian
        %   phi: constraint values
        %   cInfo: contact information
        %   settings: display and system settings ?
        %   cache: cache of quickSolve to recycle values/warm start
        %   td: timing data
    end
end

