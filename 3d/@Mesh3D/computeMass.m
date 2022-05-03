function m = computeMass(mesh)
    % computeMass computes the diagonal of the mass matrix
    %   computeMass( mesh, rho ) returns the lumped node diagonal of the mass 
    %   matrix for the given mesh

    m = zeros(mesh.N * 3, 1);
    el = mesh.el;
    for i = 1:size(el, 1)
        rho = mesh.materials(mesh.materialIndex(i)).rho;
        for j = 1:4
            %mass is the density multiplied by the area divided by the
            %number of nodes of the element
            m(el(i).t(j) * 3 - 2) = m(el(i).t(j) * 3 - 2) + 0.25 * el(i).area * rho;
            m(el(i).t(j) * 3 - 1) = m(el(i).t(j) * 3 - 1) + 0.25 * el(i).area * rho;
            m(el(i).t(j) * 3)     = m(el(i).t(j) * 3)     + 0.25 * el(i).area * rho;
        end
    end
end