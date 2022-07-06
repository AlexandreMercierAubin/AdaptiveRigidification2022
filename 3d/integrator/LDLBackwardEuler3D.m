classdef LDLBackwardEuler3D < Integrator
    % CHOLESKIBACKWARDEULER Backward euler integrator that uses choleski to
    % solve for M - h^2K

    properties
        % none needed
    end
    
    methods
        function obj = LDLBackwardEuler3D()
            obj@Integrator();
            obj.Name = 'Backward Euler for 3D';
        end
        
        function integrate( obj, mesh, h, Jc, phi, cInfo, settings, cache, td, animationScripter, frame)
            if nargin < 4
                Jc = zeros( 0, mesh.N*3 ); % no constraints
            end          
            if nargin < 5
                phi = [];
            end
            
            ticIntegrateForces = tic;
            
            [bigB, bigGamma] = mesh.getB(cache);
            Jc = Jc * bigGamma;
            
            %F = mesh.B( mesh.ActiveBRows, : ) * mesh.p;
            F = cache.F( mesh.ActiveBRows );
            
            mesh.resetForce;
            
             % Gravity (z is vertical in matlab)
            mesh.applyAcceleration([0; 0; obj.Gravity]);
            mesh.f = mesh.f + cache.elasticForces + obj.CustomForce;
            
            %scripted animations force impulses
            for i=1:numel(animationScripter)
                mesh.f = animationScripter{i}.scriptForceAnimation(mesh.f, frame, h);
            end
            
            if isa( mesh, 'AdaptiveMesh3D' )
                mesh.computeRigidForce();
            end
            
            % compute right hand side
            rhs = h * mesh.getCurrentForce;

            C = cache.C( mesh.ActiveBRows, mesh.ActiveBRows );
            bigAlpha1 = mesh.bigAlpha1( mesh.ActiveBRows, mesh.ActiveBRows );
            
            % assembly is somewhat costly... perhaps can be a bit more
            % efficient about it...
            % 1) sparsity structure is fixed... this can be exploited
            % 2) K and Kd do not need to be built separately!
            
            K = sparse(bigB' * C * bigB); % sparse here is in case there's no C (because all rigid)
            Kd = sparse(bigB' * (bigAlpha1 * C) * bigB); % sparse here is in case there's no C (because all rigid)
            
            v = mesh.getCurrentVelocity();
            Md = bigGamma' * mesh.Md * bigGamma;
            M = mesh.getM();
      
            rhs = rhs - h*Md*v + h*Kd*v + h^2*K*v;
            A = M - h * (-Md + Kd) - h^2 * K;
            
            % make set of total active dofs
            ii = mesh.activeDOFs;
            [L, D, P, S] = ldl(A(ii, ii));
                        
            % part of deltav that comes from outside forces
            % this is basically Ainv * rhs but without computing Ainv
            extDv = S * (P * (L' \ (D \ (L \ (P' * (S * rhs(ii)))))));

            td.integrateForces = toc( ticIntegrateForces );
            ticIntegrateContact = tic;
            if isa(mesh, "AdaptiveMesh3D")
                deltav = zeros( 3 * size(mesh.ElasticInds, 2) + 6*numel(mesh.RigidBodies), 1);
            else
                deltav = zeros( mesh.N*3, 1);
            end
            
            if ( isempty(cInfo) || numel(ii)==0 )
                deltav(ii) = extDv;
            else            
                Lrhs = Jc(:, ii) * (v(ii)+extDv); % lower right hand side           
                %adding baumgarte to normal velocities
                Lrhs(1:3:end) = Lrhs(1:3:end) + obj.Baumgarte * phi;
                %but not the the friction
                % LRHS gets velocity term for contacts with moving obstacles
                for i = 1:numel(cInfo)
                    Lrhs(i*3-2:i*3) = Lrhs(i*3-2:i*3) + cInfo(i).velocity;
                end

                n = size(Lrhs, 1);
                warmStartLambdas = zeros(n,1);
                if ( settings.WarmStartEnabled )
                    if ~isempty( cache.prevCInfo )
                        warmStartLambdas = cache.findWarmStartLambdaArray();
                    end
                end
 

                assert( issparse( Jc ) );
%                 [lambda, dv] = solveLDLTPGS3D(settings.PGSiterations, Jc(:, ii), L, D, P, S, Lrhs, warmStartLambdas, cInfo, obj.Compliance, td);
                % The with JAinvJT version is more than 2x faster... not
                % clear it ever makes sense to use the T version commented
                % above.
                [lambda, dv] = solveLDLTPGS3DwithJAinvJT(settings.PGSiterations, Jc(:, ii), L, D, P, S, Lrhs, warmStartLambdas, cInfo, obj.Compliance, td);

                cache.prevLambdas = lambda;

                % final deltav is solved deltav from contacts + deltav from
                % external forces
                deltav(ii) = dv + extDv;
            end

            % store contact information for warm starts on next run
            cache.prevCInfo = cInfo;
            % the warm start cache has been used, can clear it now
            cache.clearWarmStartInfo;
            
            td.integrateContacts = toc( ticIntegrateContact );

            cache.oldp = mesh.p;
            cache.oldv = mesh.v;
            
            mesh.updateParticles(h, deltav(ii));
        end
    end
end