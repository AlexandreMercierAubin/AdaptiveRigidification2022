classdef BackwardEuler3D < Integrator
    % CHOLESKIBACKWARDEULER Backward euler integrator that uses choleski to
    % solve for M - h^2K

    properties
        % none needed
        prevWarmstarts
        regularizator = 1;
        useGamma = 0;
    end
    
    methods
        function obj = BackwardEuler3D()
            obj@Integrator();
            obj.Name = 'Backward Euler for 3D';
        end
        
        function integrate( obj, mesh, h, Jc, phi, cInfo, settings, cache, td, animationScripts, frame  )
            if nargin < 4
                Jc = zeros( 0, mesh.N*3 ); % no constraints
            end          
            if nargin < 5
                phi = [];
            end
            
            ticIntegrateForces = tic;
            ii = mesh.activeDOFs;
            
            [bigB, bigGamma] = mesh.getB(cache);
            JcAdaptive = sparse(Jc * bigGamma);
            
            mesh.resetForce;
            
             % Gravity (z is vertical in matlab)
            mesh.applyAcceleration([0; 0; obj.Gravity]);
            mesh.f = mesh.f + cache.elasticForces + obj.CustomForce;
            
            for i=1:numel(animationScripts)
                mesh.f = animationScripts{i}.scriptForceAnimation(mesh.f, frame, h);
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
            
%             Ainv = inv(A(ii,ii));
            Aii = A(ii,ii);
%             extDeltav = Ainv *rhs(ii);
            extDeltav = Aii \ rhs(ii);
                        
            % part of deltav that comes from outside forces
            % this is basically Ainv * rhs but without computing Ainv
            
            td.integrateForces = toc( ticIntegrateForces );
            ticIntegrateContact = tic;
            if isa(mesh, "AdaptiveMesh3D")
                deltav = zeros( 3 * size(mesh.ElasticInds, 2) + 6*numel(mesh.RigidBodies), 1);
            else
                deltav = zeros( mesh.N*3, 1);
            end
            
            if ( isempty(cInfo) )
                deltav(ii) = extDeltav;
            else            
                Lrhs = JcAdaptive(:, ii) * (v(ii) + extDeltav); % lower right hand side         
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
                        [warmStartLambdas,newContacts] = cache.findWarmStartLambdaArray();
                    end
                end
 
%                 Acontact = -JcAdaptive(:, ii)*Ainv*JcAdaptive(:, ii)';
                Acontact = -JcAdaptive(:, ii)*(Aii\JcAdaptive(:, ii)');
                if obj.regularizator == 1
                    compliance = eye(size(Jc,2));
                    Winv = Jc*compliance*Jc';
                elseif obj.regularizator == 2
                    W = mesh.Mii;
                    Winv = Jc*(W\Jc');
                elseif obj.regularizator == 3
                    W = cache.Lpre*cache.Lpre';
                    Winv = Jc*(W\Jc');
                elseif obj.regularizator == 4
                    W = cache.AInvBlocks;
                    Winv = Jc*(W*Jc');
                elseif obj.regularizator == 5
                    W = cache.Apre;
                    Winv = Jc*(W\Jc');
                elseif obj.regularizator == 6
                    elC = cache.C;
                    elAlpha1 = mesh.bigAlpha1;
                    elK = sparse(mesh.B' * elC * mesh.B);
                    elKd = sparse(mesh.B' * (elAlpha1 * elC) * mesh.B);
                    elM = mesh.M;

                    dampingTerm = h*(-mesh.Md + elKd);
                    W = elM  - h^2 * elK - dampingTerm;
                    Winv = Jc*(W\Jc');
                end
                
                if obj.useGamma == 1
                    gamma = obj.Compliance/eigs(Winv,1);
                elseif obj.useGamma == 0
                    gamma = 1;
                else
                    gamma = obj.useGamma;
                end

                Acontact = Acontact - gamma*Winv;
                
                obj.prevWarmstarts = warmStartLambdas;
                
                assert( issparse( Jc ) );
                lambda = solvePGS3D(settings.PGSiterations, Acontact, Lrhs, warmStartLambdas, cache.cInfo);

                cache.prevLambdas = lambda;

                % final deltav is solved deltav from contacts + deltav from
                % external forces
                cache.prevLambdas = lambda;
%                 deltav(ii) = Ainv *(JcAdaptive(:, ii)' * lambda) + extDeltav;
                deltav(ii) = (Aii \ (JcAdaptive(:, ii)' * lambda)) + extDeltav;
            end

            % store contact information for warm starts on next run
            cache.prevCInfo = cInfo;
            % the warm start cache has been used, can clear it now
            cache.clearWarmStartInfo;
            
            td.integrateContacts = toc( ticIntegrateContact );
            tic
            
            cache.oldp = mesh.p;
            cache.oldv = mesh.v;
            
%             for i = 1: numel(animationScripts)
%                 deltav = animationScripts{i}.scriptVelocity(mesh.p, mesh.v, deltav, mesh.ElasticDOFs, mesh.ActiveDofsCorrespondingID, frame, mesh.rigidIDbyVert, h);
%             end
            deltav = deltav(ii);
            
            mesh.updateParticles(h, deltav);
            cache.oldDv = mesh.v - cache.oldv;
            update = toc;
            td.integrateForces = td.integrateForces + update;
        end
    end
end