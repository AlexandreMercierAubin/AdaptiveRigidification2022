function [comparisons, initialMeshes, meshes, integrators, rigidificators] = setupComparisons(meshes, integrators, rigidificators)
   %[comparisons, initialMeshes, meshes, integrators, rigidificators] = setupComparisons(meshes, integrators, rigidificators)
   %Make the scene comparisons so we can overlap meshes with different
   %scene settings for instance materials or adaptive rigidification
   %thresholds.

    % number of simulations to compare
    comparisons = max([numel(meshes), numel(integrators), numel(rigidificators)]);
    
    % detect invalid configurations
    invalidMeshes = numel(meshes) ~= 1 && numel(meshes) ~= comparisons;
    invalidIntegrators = numel(integrators) ~= 1 && numel(integrators) ~= comparisons;
    invalidRigidificators = numel(rigidificators) ~= 1 && numel(rigidificators) ~= comparisons;
    %invalidContacts = numel(contactFinders) ~= 1 && numel(contactFinders) ~= comparisons;
    
    if invalidMeshes || invalidIntegrators || invalidRigidificators
        disp('Error: mismatch between number of mesh sets, integrators, rigidificators and contact finders');
        return;
    end
    
    % copy missing parameters to make them all match the number of comparisons
    if numel(meshes) == 1
        for k = 2:comparisons
            meshes{k} = meshes{1}.clone();
        end
    end
    
    if numel(integrators) == 1
        for k = 2:comparisons
            integrators{k} = integrators{1};
        end
    end
    
    if numel(rigidificators) == 1
        for k = 2:comparisons
            rigidificators{k} = rigidificators{1};
        end
    end

    initialMeshes = cell(numel(meshes),1);
    for i = 1:comparisons
        initialMeshes{i}=meshes{i}.clone();
    end
end

