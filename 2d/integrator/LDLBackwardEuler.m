classdef LDLBackwardEuler < Integrator
    % CHOLESKIBACKWARDEULER Backward euler integrator that uses choleski to
    % solve for M - h^2K

    properties
        spyAfig
        ldlType = 3;
    end
    
    methods
        function obj = LDLBackwardEuler()
            obj@Integrator();
            obj.Name = 'Choleski Backward Euler';
        end
        
        function integrate( obj, mesh, h, Jc, phi, cInfo, settings, cache, td, animationScripts, frame )
            
            if nargin < 4
                Jc = sparse(zeros( 0, mesh.N*2 )); % no constraints
            end          
            if nargin < 5
                phi = [];
            end
            
            ticIntegrateForces = tic;
            
            [bigB, bigGamma] = mesh.getB(cache);

            % convert Jc from all dofs to elastic + rigid dofs
            Jac = sparse(Jc * bigGamma);
            
            mesh.resetForce;
            
            mesh.applyAcceleration([0; obj.Gravity]);

            % Here is how the paper does it, using only active B rows and
            % with only the dpsidF for those elastic elements (sigma in the
            % paper). I can't see the difference.  Compute time might not
            % be saved by doing this as we alrady have the elastic forces
            % for the quicksolve.
%              BE = mesh.B(mesh.ActiveBRows, :); % B but only elastic element rows 
%              elasticForcesNonRigidOnly = BE' * cache.dpsidF( mesh.ActiveBRows ); 
%              mesh.f = mesh.f + elasticForcesNonRigidOnly + obj.CustomForce;

            mesh.f = mesh.f + cache.elasticForces + obj.CustomForce;
            
            %add animationScript forces when scripted
            for i=1:numel(animationScripts)
                mesh.f = animationScripts{i}.scriptForceAnimation(mesh.f, frame, h);
            end

            % see function below, computes the forces on each rigid bodies
            % given the force on each particle
            if isa( mesh, 'AdaptiveMesh' )
                mesh.computeRigidForce();
            end
            %% big force column vector being put in the RHS
            % get current force is different than mesh.f as it will include
            % linear and rotational rigid body forces for adaptive meshes
            rhsAdaptive = h * mesh.getCurrentForce; 
            
            Md = mesh.Md;
            
            bigCAdaptive = cache.C( mesh.ActiveBRows, mesh.ActiveBRows );
            bigAlpha1G = mesh.bigAlpha1( mesh.ActiveBRows, mesh.ActiveBRows );
            
            KAdaptive = sparse(bigB' * bigCAdaptive * bigB); % sparse here is in case there's no C (because all rigid)
            KdAdaptive = sparse(bigB' * (bigAlpha1G * bigCAdaptive) * bigB); % sparse here is in case there's no C (because all rigid)
            
            % note that getM returns the adaptive mass matrix if the mesh
            % is adaptive, meaning that it includes the mass and moment of
            % inertia of the rigid bodies
            MAdaptive = mesh.getM();           
            
            % build big velocity column vector and add rayleigh damping to
            % the RHS. Note tha tin calling getCurrent velocity it isn't
            % must mesh.v but will contain the rigid velocities too.
            vG = mesh.getCurrentVelocity();
            
            % mass damping matrix Md 
            MdAdaptive = bigGamma' * Md * bigGamma;
                        
            rhsAdaptive = rhsAdaptive - h*MdAdaptive*vG + h*KdAdaptive*vG + h^2*KAdaptive*vG;
            
            AAdaptive = MAdaptive - h * (-MdAdaptive + KdAdaptive) - h^2 * KAdaptive; 

            % LDL avoids building A^-1
            ii = mesh.activeDOFs;
      
            % part of deltav that comes from outside forces
            % this is basically Ainv * rhs but without computing Ainv
            if obj.ldlType == 1 %The types are mostly there to test the efficiency of the variations of the decompositions
                [LG, DG] = ldl(full(AAdaptive(ii, ii)));  
                a = LG\rhsAdaptive(ii);
                b = DG\a;
                extDvG = LG'\b;
            elseif obj.ldlType == 2
                [LG, DG, PG] = ldl(AAdaptive(ii, ii));
                extDvG = (PG * (LG' \ (DG \ (LG \ (PG' * (rhsAdaptive(ii)))))));
            else %This is the standard and efficient way of doing a LDL decomposition
                [LG, DG, PG, SG] = ldl(AAdaptive(ii, ii));
                extDvG = SG * (PG * (LG' \ (DG \ (LG \ (PG' * (SG * rhsAdaptive(ii)))))));
            end

            td.integrateForces = toc( ticIntegrateForces );

            ticIntegrateContact = tic;

            if isa(mesh,"AdaptiveMesh")
                deltav = zeros( 2*numel(mesh.ElasticInds) + numel(mesh.RigidBodies) * 3, 1);
            else
                deltav = zeros( mesh.N*2, 1);
            end
            
            if ( isempty(cInfo) )
                deltav(ii) = extDvG;
            else            
                LrhsAdaptive = Jac(:, ii) * (vG(ii) + extDvG); % lower right hand side
                %adding baumgarte to normal velocities
                LrhsAdaptive(1:2:end) = LrhsAdaptive(1:2:end) + obj.Baumgarte * phi;
                %but not the the friction
                % LRHS gets velocity term for contacts with moving obstacles
                for i = 1:numel(cInfo)
                    LrhsAdaptive(i*2-1:i*2) = LrhsAdaptive(i*2-1:i*2) + cInfo(i).velocity;
                end
                
                n = size(LrhsAdaptive, 1);
                warmStartLambdas = zeros(n,1);
                if ( settings.WarmStartEnabled )
                    if ~isempty( cache.prevCInfo )
                        warmStartLambdas = cache.findWarmStartLambdaArray();
                    end
                end
 
                assert( issparse( Jac ) );
               
                % The withJAinvJT version is faster for the case when there are
                % tons of DOFs compared to the number of contacts, and
                % it seems that this is also true for the case when
                % there is only a single rigid body too.
                % (about 25% faster for ground_iShapeHD... for
                % instance)

                [lambda, dv] = solveLDLTPGSwithJAinvJT(settings.PGSiterations, Jac(:, ii), LG, DG, PG, SG, LrhsAdaptive, warmStartLambdas, cInfo, obj.Compliance );

%                 This residual is missing the compliance term
                %residual = sum(Jac(:,ii)*AAdaptive(ii,ii)*Jac(:,ii)'*lambda - LrhsAdaptive)
                
                cache.prevLambdas = lambda;

                % final deltav is solved deltav from contacts + deltav from
                % external forces
                deltav(ii) = dv + extDvG;
            end
            
            
            % store contact information for warm starts on next run
            cache.prevCInfo = cInfo;
            % the warm start cache has been used, can clear it now
            cache.clearWarmStartInfo;
            
            td.integrateContacts = toc( ticIntegrateContact );
            
            tic
            
            cache.oldp = mesh.p;
            cache.oldv = mesh.v;
            
            deltav = deltav(ii);
            
            mesh.updateParticles(h, deltav);
            update = toc;
            td.integrateForces = td.integrateForces + update;

        end
    end
end