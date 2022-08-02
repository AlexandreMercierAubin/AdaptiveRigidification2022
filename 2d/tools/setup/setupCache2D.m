function [caches] = setupCache2D(settings, meshes, rigidificators, comparisons, energyModel, h)
    % setupCache2D(settings, meshes, rigidificators, comparisons, h)
    % create the cache for the scene. This includes the precomputation of
    % the A inverse blocks for the quicksolve

    sceneDescriptor = [string(settings.SceneName)+string(energyModel.name)+strrep(num2str(h),'.','')];
    cachedAinvString = ['2d/data/cached/',convertStringsToChars(sceneDescriptor),'Ainv.mat'];
    toRecompute = ~isfile(cachedAinvString) || settings.recomputeCacheAinv;
    if ~toRecompute
        cachestmp = load(cachedAinvString);
        if (numel(cachestmp.caches) ~= comparisons) && ~settings.sameComparisonAinv
            toRecompute = true;
        else
            caches = cachestmp.caches;
            if settings.sameComparisonAinv && numel(caches) < comparisons
                for k = numel(caches):comparisons
                    caches{k} = caches{1};
                    copyStream = getByteStreamFromArray(caches{1});
                    caches{k} = getArrayFromByteStream(copyStream);
                 end
                save(cachedAinvString,'caches');
            end
        end
    end
    
    if toRecompute
        caches = cell(comparisons, 1);
        for k = 1:comparisons
            caches{k} = QuickSolveCache();
            caches{k}.ActiveB = meshes{k}.B;
            if settings.sameComparisonAinv && k > 1
                caches{k} = caches{1};
                copyStream = getByteStreamFromArray(caches{1});
                caches{k} = getArrayFromByteStream(copyStream);
                continue;
            end
            
            rigidificator = rigidificators{k};
            mesh2D = meshes{k};
            if ( ~isa(mesh2D, 'AdaptiveMesh') )
                continue;
            end
            cache = caches{k};
            cache.oldv = mesh2D.v;
            cache.oldp = mesh2D.p;
            precomputeAInverseDiagonalBlocks( cache, mesh2D, h , energyModel);

            %preconditioner for quicksolve
            %compute young's moduli
            mu = mesh2D.elMu;
            lam = mesh2D.elLambda;
            %eq from https://subsurfwiki.org/wiki/Young's_modulus in other
            %expressions
            youngs_moduli = (mu.*(3.*lam+2.*mu))./(lam+mu);
            %get to the bottom of this mystery!
            diagC = spdiags(reshape(...
                            repmat(youngs_moduli, 4,1),...
                                    4*numel(youngs_moduli),1 ...
                                ),...
                                [0],... 
                                4*numel(youngs_moduli),...
                                4*numel(youngs_moduli)...
                            );
            Kapprox = -mesh2D.B'*diagC * mesh2D.B;
            
            Kd = sparse(mesh2D.B' * (mesh2D.bigAlpha1 * diagC) * mesh2D.B);
%             dampingTerm = h * (-mesh2D.Md + Kd);
            Apre = mesh2D.M - h*h*Kapprox; %+ dampingTerm;

            %pinned vertices
            ii = mesh2D.unpinnedDOFs;

            % prep the matrix and permutation
            Aii = Apre(ii,ii);
            I = speye(size(Aii));
            if strcmpi(rigidificator.Permutation, 'SYMAMD')
                perm = symamd(Aii);
            elseif strcmpi(rigidificator.Permutation, 'SYMRCM')
                perm = symrcm(Aii);
            elseif strcmpi(rigidificator.Permutation, 'COLPERM')
                perm = colperm(Aii);
            elseif strcmpi(rigidificator.Permutation, 'COLAMD')
                perm = colamd(Aii);
            elseif strcmpi(rigidificator.Permutation, 'RANDPERM')
                perm = randperm(size(Aii,1));
            elseif strcmpi(rigidificator.Permutation, 'DISSECT')
                perm = dissect(Aii);
            elseif strcmpi(rigidificator.Permutation, 'AMD')
                perm = amd(Aii);
            else % strcmpi(rigidificator.Permutation, 'NONE')
                perm = 1:size(Aii,1);
            end
            Perm = I(perm,:);
            App = Aii(perm,perm);
            cache.App = App;
            cache.Apre = Apre;  % NOTE: this is only used for regularization testing (and is perhaps used wrong incorrectly?)
            
            % note: default diagcomp is zero for ichol, and should probably keep it that way !
            % that is, never use diagcomp, it tends to make the quicksolve pcg work poorly
            if strcmpi(rigidificator.Preconditionner, 'ICHOLICT')
                Lpre = ichol( App, struct('type','ict','droptol',1e-6,'michol','on','diagcomp', 0.0) ); 
                fprintf('ICHOLICT with droptol 1e-6 has factor %g more non-zeros than IC(0)\n', ...
                    nnz(Lpre) / nnz( ichol( App, struct('type','nofill','michol','on') ) ) );
            elseif strcmpi(rigidificator.Preconditionner, 'ICHOL')
                Lpre = ichol( App, struct('type','nofill','michol','on') ); 
            elseif strcmpi(rigidificator.Preconditionner, 'CHOL')
                Lpre = chol( App, 'lower' );
            elseif strcmpi(rigidificator.Preconditionner, 'LDL')
                Lpre = chol( App, 'lower' );   
            else
                Lpre = I(perm,perm);
            end
            
            cache.preconditioner = @(r) Perm'*(Lpre'\(Lpre\(Perm*r)));
            cache.Perm = Perm;
            cache.Lpre = Lpre;
        end
        save(cachedAinvString,'caches');
    end
end

