function quickSolveNoConstraint( cache, integrator, mesh2D, h, Jc, phi, settings, animationScripter, frame)
    % QUICKSOLVENOCONSTRAINT Performs the single iteration of
    % preconditionned conjugate gradient required for derigidification.

    if nargin < 5
         Jc = zeros( 0, mesh2D.N*2 ); % no constraints
    end  
    ii = mesh2D.unpinnedDOFs;
    
    % reset forces
    mesh2D.f = zeros(size(mesh2D.f));
    
%     mesh.applyAcceleration([0; integrator.Gravity]);
%      lets use forward euler on g to separate the gravity forces
    if integrator.separateQuicksolveGravity && cache.recomputeGravityDv == 1
        %TODO: This needs to be updated on mass change or when dofs are pinned
        %TODO: separating gravity from the pcg solve does not work with
        %pinned dofs; fix it
        gravity = zeros( mesh2D.N*2, 1 );
        gravity(2:2:end) = integrator.Gravity;
        
        cache.gravityDv = h*gravity(ii);
        cache.recomputeGravityDv = 0;
    elseif ~integrator.separateQuicksolveGravity
        mesh2D.f(2:2:end) = mesh2D.mass(2:2:end) * integrator.Gravity;
        cache.gravityDv = zeros(size(mesh2D.f(ii)));
    end
    
    mesh2D.f = mesh2D.f + cache.elasticForces + integrator.CustomForce;
    
    %add animationScript forces when scripted
    for i=1:numel(animationScripter)
        mesh2D.f = animationScripter{i}.scriptForceAnimation(mesh2D.f, frame, h);
    end
    
    rhs = h * mesh2D.f(ii);
    
    bigB = mesh2D.Bii;
    bigC = cache.C; %computeC( Fs, mesh, 1:size(mesh.t, 1) );
    bigAlpha1 = mesh2D.bigAlpha1;
    M = mesh2D.Mii;
    Md = mesh2D.Mdii;
    v = mesh2D.v(ii) + cache.gravityDv;
    
    % DONT ASSEMBLE... this will be faster
%     rhs = rhs - h*Md*v + h*Kd*v + h^2*K*v;
    % ... and while it might not be pretty, it is likewise faster to reuse
    % the common CBv expression...
    %note that this is redundant to the first call to qspcghelper within
    %the pcg
    CBv = bigC*(bigB*v);
    alphaCBv = h*(bigAlpha1*CBv) + h^2*CBv; % hKdv + h^2Kv
    fext = h*(-Md*v) + bigB'*alphaCBv; % h*(-Md*v + hKdv) +  + h^2K
    rhs = rhs + fext;
    
    % If we don't assemble K, then a multiplication by K might be
    % simpler!  Probably unavoidable to recompute all the C, but they at
    % least can be reused (those that are elastic at least)
    
    if ~integrator.useFullAinv
        Afun = @(v) qspcgHelper( h, M, Md, bigB, bigAlpha1, bigC, v);
    else
        K  = sparse(bigB' * bigC * bigB);
        Kd = sparse(bigB' * (bigAlpha1 * bigC) * bigB);
        A = M - h * (-Md + Kd) - h^2 * K;
        Afun = @(v) A*v;
    end
    
    if settings.PlotSpyA
        a = gcf;
        K  = sparse(bigB' * bigC * bigB);
        Kd = sparse(bigB' * (bigAlpha1 * bigC) * bigB);
        A = M - h * (-Md + Kd) - h^2 * K;     
        if isempty(settings.spyAfig)
            settings.spyAfig = figure;
        end
        set(0, 'CurrentFigure', settings.spyAfig);
        spy(A);
        title("spy A");
        set(0, 'CurrentFigure', a);
    end
    
    dofCount = 2 * sum([mesh2D.N]);
    
    cache.ApproximatedDeltaV = zeros(dofCount, 1);
        
    % presolve on new contacts if we have any!
    if ~integrator.PGSquicksolve
        if ~isempty( cache.cInfo ) 

%         We'll do the warm start here, and reuse the result later! This 
%         will let us know what is new, and what forces to put on the rhs
%                 
        [warmStartLambdas, isNewContact] =  cache.findWarmStartLambdaArray;
%         
%         Hmm... (sign may not make complete sense, but probably correct
%         based on what was done previously... note that we only use deltav
%         as the result of the PGS solve, so this is the main place where
%         lambdas are needed)
        tmpRhs = Jc' * warmStartLambdas;
        rhs = rhs + tmpRhs(ii);  % existing forces on the RHS

            if  any( isNewContact ) %|| true

    %             [  A  Jcn' ] [   dv   ]  = [    rhs   ]
    %             [ Jcn  0   ] [ lambda ]    [  -Jcn v  ]

    %             The key here is that we want Jcn v = 0, and solving for dv
    %             this means Jcn dv = -Jcn v for the current velocity v.
    %             So when we form the schur complement of this system we get
    %             Jcn Ainv Jcn' lambda = Jcn Ainv rhs + Jcn v 
    %             So that's what we need on the new contact rhs, ncRHS

                ncInds = find(isNewContact);
%                 ncInds = 1:size(Jc,1)/2;

    %             Create a smaller Jc matrix for new contacts
    %             Perhaps no need to solve for tangents?  only normals here!
                contactRows = ncInds*2-1; %without friction
% %                 contactRows = 1:size(Jc,1); %with friction
                Jcn = Jc( contactRows, ii );


                if ~integrator.useFullAinv
    %                 We will use precomptued diagonal 2x2 blocks of A rather than 
    %                 the proper A inverse (or LDLT factorization)
                    Ainvii = cache.AInvBlocks(ii,ii);
                    JcnAinvii = Jcn * Ainvii;
                    AtoSolve = JcnAinvii * Jcn';
                    ncRHS = JcnAinvii * rhs + Jcn * v;  
                else
    %                  -----debugging version note that this it can be more
    %                 accurate, yet slower
                    JcnAinv = Jcn/A;
                    ncRHS = JcnAinv*rhs + Jcn*v;
                    AtoSolve = JcnAinv*Jcn';
    %                 -----debugging version end
                end

    %             small and dense (only new contacts) so just do a direct solve
    %             or just let Matlab decide what's best with mldivide.
    %             newLambdas = solvePGS( 1000, AtoSolve, ncRHS, zeros(size(ncRHS)), cache.cInfo );
                newLambdas = AtoSolve \ ncRHS; 

                if (settings.DrawLambdas && settings.quicksolveSimulation)%|| settings.warmStartFromQuickSolve
                    warmStartLambdas(contactRows) = newLambdas;
                    cache.prevLambdas = warmStartLambdas;
                end
    %             NOTE: h is baked into lambda!  See the equation above and 
    %             note we are solving for dv, NOT acceleration, so lambda is
    %             an impulse!
                rhs = rhs - Jcn' * newLambdas;
            end
        end
    end
    % final deltav is solved deltav from contacts + deltav from external forces.
    
    % THIS MUST ALWAYS BE PCG IF WE WANT INFORMATION TO TRAVEL ACROSS
    % THE MESH WITH ONLY ONE ITERATION
    if settings.PCGiterations == 1 
        cache.ApproximatedDeltaV(ii) = PCGnoLoop(Afun, rhs, cache, mesh2D.objectDOFsPinned);
%       [cache.ApproximatedDeltaV(ii),~] = pcg( Afun, rhs,1e-8, settings.PCGiterations, cache.preconditioner);
    else %This else should only be called for debugging purpose
%         cache.ApproximatedDeltaV(ii) = pcg( Afun, rhs,1e-8, settings.PCGiterations, cache.preconditioner);
        cache.ApproximatedDeltaV(ii) = solvePCG( settings.PCGiterations, Afun, rhs, zeros(size(rhs)), cache );
    end
%     cache.ApproximatedDeltaV(ii) = A\rhs;

    if integrator.PGSquicksolve
        quickSolveNoConstraint();
    end 
    
    %adding the gravity if separated
     cache.ApproximatedDeltaV(ii) = cache.ApproximatedDeltaV(ii) + cache.gravityDv;
    
    for i=1:numel(animationScripter)
        cache.ApproximatedDeltaV = animationScripter{i}.scriptVelocityQuicksolve(mesh2D.p, mesh2D.v, cache.ApproximatedDeltaV, frame, h);
    end

    if settings.quicksolveSimulation
        mesh2D.v = mesh2D.v + cache.ApproximatedDeltaV;
        mesh2D.v(mesh2D.pinnedDOFs) = mesh2D.v(mesh2D.pinnedDOFs) * 0;
        mesh2D.p = mesh2D.p + h * mesh2D.v;
        cache.clearWarmStartInfo;
    end
%     semilogy(0:numel(rv)-1,rv/norm(rhs(ii)),'-o');
    %figure(2);plot(cache.DeltaV);
end