function [caches] = setupCache(settings, meshes, rigidificators, comparisons, energyModel, h)
% [caches] = setupCache(settings, meshes, rigidificators, comparisons, energyModel, h)
% create the cache for the scene. This includes the precomputation of
% the A inverse blocks for the quicksolve

sceneDescriptor = [string(settings.SceneName)+string(energyModel.name)+strrep(num2str(h),'.','')];
cachedAinvString = ['3d/data/cached/', convertStringsToChars(sceneDescriptor), 'Ainv.mat'];
    toRecompute = ~isfile(cachedAinvString) || settings.recomputeCacheAinv;
    if ~toRecompute
        cachestmp = load(cachedAinvString);
        if numel(cachestmp.caches) ~= comparisons && ~settings.sameComparisonAinv
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
            caches{k} = QuickSolveCache3D();
            caches{k}.ActiveB = meshes{k}.B;
            if settings.sameComparisonAinv && k > 1
                caches{k} = caches{1};
                copyStream = getByteStreamFromArray(caches{1});
                caches{k} = getArrayFromByteStream(copyStream);
                continue;
            end
            mesh3D = meshes{k};
            rigidificator = rigidificators{k};
            cache = caches{k};
            if ( isa(mesh3D, 'AdaptiveMesh3D') )   
                precomputeAInverseDiagonalBlocks3D( cache, mesh3D, h, energyModel);

                cache.oldv = mesh3D.v;
                cache.oldp = mesh3D.p;

                %preconditioner for quicksolve
                %compute young's moduli
                G = mesh3D.elMu;
                lam = mesh3D.elLambda;
                youngs_moduli = (G.*(6.*lam+2.*G))./(lam+G);
                Apre = mesh3D.M + h.*h.*...
                             mesh3D.B'*...
                             spdiags(reshape(repmat(youngs_moduli, 9,1), 9*numel(youngs_moduli),1), [0], 9*numel(youngs_moduli), 9*numel(youngs_moduli))*...
                             mesh3D.B;
                cache.Apre = Apre;
                
                %pinned vertices
                ii = mesh3D.unpinnedDOFs;

                % prep the matrix and permutation
                Aii = Apre(ii,ii);
                
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

                App = sparse(Aii(perm,perm));
                I = speye(size(Aii));
                Perm = I(perm,:);
                cache.App = App;
                
                setupPre = tic;
                % note: default diagcomp is zero for ichol, and should probably keep it that way !
                % that is, never use diagcomp, it tends to make the quicksolve pcg work poorly
                if strcmpi(rigidificator.Preconditionner, 'ICHOLICT')
                    Lpre = ichol( App, struct('type','ict','droptol',1e-6,'michol','on','diagcomp', 0.0) ); 
%                     fprintf('ICHOLICT with droptol 1e-6 has factor %g more non-zeros than IC(0)\n', ...
%                         nnz(Lpre) / nnz( ichol( App,
%                         struct('type','nofill','michol','on') ) ) ); % Go
                elseif strcmpi(rigidificator.Preconditionner, 'ICHOL')
                    Lpre = ichol( App, struct('type','nofill','michol','on') ); 
                elseif strcmpi(rigidificator.Preconditionner, 'CHOL')
                    Lpre = chol( App, 'lower' );
                elseif strcmpi(rigidificator.Preconditionner, 'LDL')
                    Lpre = chol( App, 'lower' );   
                else
                    Lpre = I(perm,perm);
                end
                setupPre = toc(setupPre);
                
                cache.preconditioner = @(r) Perm'*(Lpre'\(Lpre\(Perm*r)));
                cache.Perm = Perm;
                cache.Lpre = Lpre;
            end
        end
        save(cachedAinvString,'caches');
    end
end

