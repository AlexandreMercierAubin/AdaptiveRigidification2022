function quickSolveContactDebugger3D(A, Jc, integrator, phi, v,cache, settings)
    if ~isempty( cache.cInfo )
        Lrhs = Jc(:, ii) * cache.ApproximatedDeltaV(ii);
        Lrhs(1:3:end) = Lrhs(1:3:end) +  Jc(1:3:end,:) * v + integrator.Baumgarte * phi;

        Lrhs(2:3:end) = Lrhs(2:3:end) +  Jc(2:3:end,:) * v;
        Lrhs(3:3:end) = Lrhs(3:3:end) +  Jc(3:3:end,:) * v;

        for i = 1:numel(cache.cInfo)
        Lrhs(i*3-2:i*3) = Lrhs(i*3-2:i*3) + cache.cInfo(i).velocity;
        end

        n = size(Lrhs, 1);
        warmStartLambdas = zeros(n,1);
        if ( settings.WarmStartEnabled )
        if ~isempty( cache.prevCInfo )
            warmStartLambdas = cache.findWarmStartLambdaArray();
        end
        end
        %tmp test todo change for the other one
        [L, D, P, S] = ldl(A(ii, ii));
        [lambda, dv] = solveLDLTPGS3D(30, Jc(:, ii), L, D, P, S, Lrhs, warmStartLambdas, cache.cInfo, integrator.Compliance, 0);
        if settings.DrawLambdas && settings.quicksolveSimulation
         cache.prevLambdas = lambda;
        end
        %       newLambdas = solvePGS3D(1,A,rhs,warmStartLambdas,cache.cInfo);
        cache.ApproximatedDeltaV(ii) = cache.ApproximatedDeltaV(ii) + dv;
    end 
end

