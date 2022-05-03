classdef TimingData < handle
    %TIMINGDATA Keeps compute time data for comparison
    % likewise for comparison between different timing data logs
    
    properties

        % discrete quantities
        contactCount = 0;
        countParticles = 0;
        countTris = 0;
        countTotalParticles = 0;
        countTotalTris = 0;
        rigidBodies = 0;
        totalDofs = 0;
        stepSum = 0;
        labelsCountData = [ ...
                "contactCount", ...
                "countParticles", ...
                "countElements", ...
                "countTotalParticles", ...
                "countTotalElements", ...
                "rigidBodies", ...
                "totalDofs", ...
                "stepSum",...
                ];

        % time quantities
        lastContact = 0;            % collision detection
        lastSimulate = 0;           % total for the step
        fullFCForces = 0;            % computation of C forces and F        
        lastElastification = 0;
        lastElastifyQuickSolve = 0;
        integrateForces = 0;
        integrateContacts = 0;
        lastRigidification = 0;     
        SimulationLoopFull = 0;
        contactSolve1 = 0;
        contactSolve2 = 0;
        
        
        labels = [ ...
                "Collision Detection", ...
                "Simulation Step Total", ...
                "Compute Forces and C", ...
                "Elastification QuickSolve", ...
                "Elastification (Bookkeeping)", ...
                "Integrate Forces", ...
                "Integrate Contact", ...
                "Rigidification", ...
                "contact solve 1", ...
                "contact solve 2", ...
                ];
        
        % other... was for debugging and could remove (or keep)
        linearMomentum = [];
        
        log  % a matrix of #data by #frames timing data
        logCounts  % a matrix of #discretedata by #frames count data
        
        p % a matrix of #vertices by #frames of positions p at frame i. filled only if settings.RecordFramePositionInTD is true
    end
    
    methods
        function obj = TimingData()
        end
        
        function logData( obj )
            dataCounts = [
                obj.contactCount; ...
                obj.countParticles; ...
                obj.countTris; ...
                obj.countTotalParticles; ... 
                obj.countTotalTris; ...
                obj.rigidBodies; ...
                obj.totalDofs; ...
                obj.stepSum; ...
                ];
            obj.logCounts = [ obj.logCounts, dataCounts ];

            % logData Records the last stats into a matrix
            data = [ ...
                    obj.lastContact; ...            % collision detection
                    obj.lastSimulate; ...           % total for the step
                    obj.fullFCForces; ...            % computation of C forces and F        
                    obj.lastElastification; ...
                    obj.lastElastifyQuickSolve; ...
                    obj.integrateForces; ...
                    obj.integrateContacts; ...
                    obj.lastRigidification; ...  
                    obj.contactSolve1; ...
                    obj.contactSolve2; ...
                    ];
            obj.log = [ obj.log, data ];
        end
    
        function rng = getRelevantTimingIndicies( obj, varargin ) 
            % choose indicides based on timing data exceeding a threshold
            % default threshold is 1e-4, otherwise specify an argument.
            threshold = 1e-3;
            if numel(varargin) > 0
                threshold = varargin{1};
            end
            maxvals =  max( obj.log(:,5:end), [], 2 );
            [maxvalssorted, ind] = sort(maxvals, 'descend');
            rng = ind( maxvalssorted > threshold );
        end
        
        function plot( obj, varargin )
            % PLOT Draws timing data
            %    Optional arguments:
            %       Rng         an array of indices to plot, default is those things exceeding 1e-4 s
            %       LineSpec    default is '-'
            %       Smoothing   default is false, if true then smoothing filter has width 5
            %       Log         default is false
            
            LineSpec = '-';
            rng = obj.getRelevantTimingIndicies;
            plotLog = false;
            plotWithSmoothing = true;
            
            params_to_variables = containers.Map( ...
                {'Rng', 'LineSpec', 'Smoothing', 'Log'}, ...
                {'rng', 'LineSpec', 'plotWithSmoothing', 'plotLog'} );
            v = 1;
            while v <= numel(varargin)
                param_name = varargin{v};
                if isKey( params_to_variables, param_name )
                    assert(v+1<=numel(varargin));
                    v = v+1;
                    feval( @() assignin( 'caller', params_to_variables(param_name), varargin{v} ) );
                else
                    error( "Unsupported parameter: " + varargin{v} );
                end
                v=v+1;
            end
   
            % Make nice colours so we can figure out what is going on
            HSV = ones( numel(rng) + 1, 3 );
            HSV(:,1) = linspace( 0, 1, numel(rng) + 1 );
            HSV(:,2) = linspace( 1, 0.6, numel(rng) + 1 );
            HSV(:,3) = HSV(:,3)*0.75;
            RGB = hsv2rgb( HSV );

            flog1 = obj.log;    
            if ( plotWithSmoothing ) 
                numSmoothFrames = 5;
                flog1 = [filtfilt( ones(1,numSmoothFrames)./numSmoothFrames, 1, flog1' )]';
            end
            
            for i = 1:numel(rng)
                if ( plotLog )
                    semilogy( flog1( rng(i), 10:end )', LineSpec, 'Color', RGB(i,:), 'LineWidth', 2 );
                else
                    plot( flog1( rng(i), 10:end )', LineSpec, 'Color', RGB(i,:), 'LineWidth', 2 );
                end
                hold on;
            end
            legend( obj.labels( rng ) );
        end

        function plotCompare( obj, td, varargin )
            % figure out the common interesting range and pass the same 
            % range to both plots to get constistent colours!
            
            rng = union( obj.getRelevantTimingIndicies, td.getRelevantTimingIndicies );
            plotLog = false;
            plotWithSmoothing = true;
            
            params_to_variables = containers.Map( ...
                {'Rng', 'Smoothing', 'Log'}, ...
                {'rng', 'plotWithSmoothing', 'plotLog'} );
            v = 1;
            while v <= numel(varargin)
                param_name = varargin{v};
                if isKey( params_to_variables, param_name )
                    assert(v+1<=numel(varargin));
                    v = v+1;
                    feval( @() assignin( 'caller', params_to_variables(param_name), varargin{v} ) );
                else
                    error( "Unsupported parameter: " + varargin{v} );
                end
                v=v+1;
            end
            
            obj.plot( 'LineSpec', '-', 'Rng', rng, 'Log', plotLog, 'Smoothing', plotWithSmoothing );
            td.plot( 'LineSpec', ':', 'Rng', rng, 'Log', plotLog, 'Smoothing', plotWithSmoothing );
        end
        
        function plotCompareSiggraphVersion( obj, td, h )
            % figure out the common interesting range and pass the same 

            plotLog = false;
            plotWithSmoothing = true;
            
            flog1 = obj.log;  
            numSmoothFrames = 5;
            flog1 = [filtfilt( ones(1,numSmoothFrames)./numSmoothFrames, 1, flog1' )]';
            
            flog2 = td.log;  
            numSmoothFrames = 5;
            flog2 = [filtfilt( ones(1,numSmoothFrames)./numSmoothFrames, 1, flog2' )]';
            
            hold on;
            plot([h*4:h:h*size(flog1,2)],flog1(2,4:end), 'LineWidth', 2, 'Color',[170/255, 20/255, 20/255]);
            plot([h*4:h:h*size(flog2,2)],flog2(2,4:end), 'LineWidth', 2, 'Color', [204/255, 85/255, 0/255] );
            plot([h*4:h:h*size(flog2,2)],flog2(4,4:end) + flog2(5,4:end), 'LineWidth', 2, 'Color', [20/255, 20/255, 170/255] );
        end
    end
end

