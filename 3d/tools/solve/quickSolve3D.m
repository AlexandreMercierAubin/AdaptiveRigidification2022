function quickSolve3D( cache, integrator, mesh3D, h, Jc, phi, settings, animationScripter, frame)
    % QUICKSOLVE performs the single iteration of preconditionned conjugate
    % gradient required for derigidification, does not support constraints
    if nargin < 5 || isempty(Jc)
         Jc = zeros( 0, mesh3D.N*3); % no constraints
    end
    
    cgTic = tic;
     
    % NOTE: OPTIMIZATIONS...
    % quicksolve computes F and STVK forces and C for *all* elements and
    % this can be reused in the solve later (though only small parts are
    % needed).  All of this is cheap and linear in #TRI but gets expensive.
    
    bigB = mesh3D.B;
    %F = mesh.B * mesh.p;
    
    % update forces
    mesh3D.f = zeros( mesh3D.N*3, 1 );
    % Gravity (z is vertical in matlab)
    %lets use forward euler on g to separate the gravity forces
    ii = mesh3D.unpinnedDOFs;
    if integrator.separateQuicksolveGravity && cache.recomputeGravityDv == 1
        %TODO: This needs to be updated on mass change or when dofs are pinned
        %TODO: separating gravity from the pcg solve does not work with
        %pinned dofs; fix it
        gravity = zeros( mesh3D.N*3, 1 );
        gravity(3:3:end) = integrator.Gravity;
        
        cache.gravityDv = h*gravity(ii);
        cache.recomputeGravityDv = 0;
    elseif ~integrator.separateQuicksolveGravity
        mesh3D.f(3:3:end) = mesh3D.mass(3:3:end) * integrator.Gravity;
        mesh3D.f(mesh3D.pinnedDOFs) = 0;
        cache.gravityDv = zeros(size(mesh3D.f(ii)));
    end

    mesh3D.f = mesh3D.f + cache.elasticForces + integrator.CustomForce;
    
    %add animationScript forces when scripted
    for i=1:numel(animationScripter)
        mesh3D.f = animationScripter{i}.scriptForceAnimation(mesh3D.f, frame, h);
    end
    
    rhs = h * mesh3D.f(ii);
    %params = zeros(size(mesh.t, 1), 2);
    
    bigB = mesh3D.Bii;
    bigC = cache.C; %computeC( Fs, mesh, 1:size(mesh.t, 1) );
    bigAlpha1 = mesh3D.bigAlpha1;
    if integrator.useFullAinv
        K  = sparse(bigB' * bigC * bigB);
        Kd = sparse(bigB' * (bigAlpha1 * bigC) * bigB);
    end
    M = mesh3D.Mii;
    Md = mesh3D.Mdii;
    v = mesh3D.v(ii) + cache.gravityDv;
    
    % DONT ASSEMBLE... this will be faster
    %rhs1 = rhs - h*Md*v + h*Kd*v + h^2*K*v;
    % ... and while it might not be pretty, it is likewise faster to reuse
    % the common CBv expression...
    CBv = bigC*(bigB*v);
    alphaCBv = h*(bigAlpha1*CBv) + h^2*CBv;
    fext = -h*(Md*v) + bigB'*alphaCBv;
    rhs = rhs + fext;
    
    % If we don't assemble K, then a multiplication by K might be
    % simpler!  Probably unavoidable to recompute all the C, but they at
    % least can be reused (those that are elastic at least)
    
    if integrator.useFullAinv
        % This is an option to compare/test the use of the full Ainv
        % vs the unassembled version. Both should give the same solution
        A = M - h * (-Md + Kd) - h * h * K;     
        Afun = @(v) A(ii,ii)*v;
    else
        Afun = @(v) qspcgHelper3D( h, M, Md, bigB, bigAlpha1, bigC, v);
    end
    
    dofCount = 3 * sum([mesh3D.N]);
    
    cache.ApproximatedDeltaV = zeros(dofCount, 1);
    
    if ~integrator.PGSquicksolve
        % presolve on new contacts if we have any!
        if ~isempty( cache.cInfo ) 

            % We'll do the warm start here, and reuse the result later! This 
            % will let us know what is new, and what forces to put on the rhs

            [warmStartLambdas, isNewContact] =  cache.findWarmStartLambdaArray;

            % Hmm... (sign may not make complete sense, but probably correct
            % based on what was done previously... note that we only use deltav
            % as the result of the PGS solve, so this is the main place where
            % lambdas are needed)
            tmpRhs = Jc' * warmStartLambdas;
            rhs = rhs + tmpRhs(ii); % existing forces on the RHS

            if any( isNewContact )

                % [  A  Jcn' ] [   dv   ]  = [    rhs   ]
                % [ Jcn  0   ] [ lambda ]    [  -Jcn v  ]

                % The key here is that we want Jcn v = 0, and solving for dv
                % this means Jcn dv = -Jcn v for the current velocity v.
                % So when we form the schur complement of this system we get
                % Jcn Ainv Jcn' lambda = Jcn Ainv rhs + Jcn v 
                % So that's what we need on the new contact rhs, ncRHS

                ncInds = find(isNewContact);

                % Create a smaller Jc matrix for new contacts
                % Perhaps no need to solve for tangents?  only normals here!

               Jcn = Jc( ncInds*3-2, ii );


                if integrator.useFullAinv
                     % -----debugging version
                    %uses the full A
                    JcnAinv = Jcn/A;
                    ncRHS = JcnAinv*rhs + Jcn*v;
                    AtoSolve = JcnAinv*Jcn';
                else
                    % We will use precomptued diagonal 3x3 blocks of A rather than 
                    % the proper A inverse (or LDLT factorization)
                    Ainvii = cache.AInvBlocks(ii,ii);
                    JcnAinvii = Jcn * Ainvii;
                    AtoSolve = JcnAinvii * Jcn';
                    ncRHS = JcnAinvii * rhs + Jcn * v;
                end

                % small and dense (only new contacts) so just do a direct solve
                % or just let Matlab decide what's best with mldivide.
    %             newLambdas = solvePGS3D(1,AtoSolve,ncRHS,zeros(numel(ncRHS),1),cache.cInfo);
                newLambdas = AtoSolve \ ncRHS; 

                % NOTE: h is baked into lambda!  See the equation above and 
                % note we are solving for dv, NOT acceleration, so lambda is
                % an impulse!
                rhs = rhs - Jcn' * newLambdas;
            end
        end
    end

    % final deltav is solved deltav from contacts + deltav from external forces.
    % This is always 1 iteration, the other version is also there to test
    % convergence
    if settings.PCGiterations == 1 %This version should be much faster for 1 iteration by avoiding wasteful computations
        cache.ApproximatedDeltaV(ii) = PCGnoLoop(Afun, rhs, cache, mesh3D.objectDOFsPinned);
    else
        cache.ApproximatedDeltaV(ii) = solvePCG( settings.PCGiterations, Afun, rhs, zeros(size(rhs)), cache );
    end
      
    if integrator.PGSquicksolve %PGS contact handler to used only to debug contacts in the quicksolve
        quickSolveContactDebugger3D(A, Jc, integrator, phi, v,cache, settings)
    end
    
    % update the change in velocities
    cache.ApproximatedDeltaV(ii) = cache.ApproximatedDeltaV(ii) + cache.gravityDv;

    %add animations
    for i=1:numel(animationScripter)
        cache.ApproximatedDeltaV = animationScripter{i}.scriptVelocityQuicksolve(mesh3D.p, mesh3D.v, cache.ApproximatedDeltaV, frame, h);
    end

    %debug option that allows to visualize the results of the quicksolve
    if settings.quicksolveSimulation %
        mesh3D.v = mesh3D.v + cache.ApproximatedDeltaV;
        mesh3D.v(mesh3D.pinnedDOFs) = mesh3D.v(mesh3D.pinnedDOFs) * 0;
        mesh3D.p = mesh3D.p + h * mesh3D.v;
        cache.clearWarmStartInfo();
    end
    
end