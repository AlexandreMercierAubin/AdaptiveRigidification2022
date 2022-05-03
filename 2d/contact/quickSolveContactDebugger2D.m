function quickSolveContactDebugger2D(A, Jc, integrator, phi, v,cache, settings)
    if ~isempty( cache.cInfo )
        Lrhs = Jc(:, ii) * cache.ApproximatedDeltaV(ii);
        Lrhs(1:2:end) = Lrhs(1:2:end) +  Jc(1:2:end,ii) * v + integrator.Baumgarte * phi;
        %but not the the friction
        Lrhs(2:2:end) = Lrhs(2:2:end) +  Jc(2:2:end,ii) * v;
        % LRHS gets velocity term for contacts with moving obstacles
        for i = 1:numel(cache.cInfo)
            Lrhs(i*2-1:i*2) = Lrhs(i*2-1:i*2) + cache.cInfo(i).velocity;
        end

        n = size(Lrhs, 1);
        warmStartLambdas = zeros(n,1);
        if ( settings.WarmStartEnabled )
            if ~isempty( cache.prevCInfo )
                warmStartLambdas = cache.findWarmStartLambdaArray();
            end
        end

        [L, D, P, S] = ldl(A);
        [lambda, dv] = solveLDLTPGS(100, Jc(:, ii), L, D, P, S, Lrhs, warmStartLambdas, cache.cInfo, integrator.Compliance, 0);

        if settings.DrawLambdas && settings.quicksolveSimulation
             cache.prevLambdas = lambda;
        end
        % final deltav is solved deltav from contacts + deltav from
        % external forces
        cache.ApproximatedDeltaV(ii) = cache.ApproximatedDeltaV(ii) + dv;
    end 
end

