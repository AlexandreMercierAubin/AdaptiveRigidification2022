function m = computeMass(mesh2d)
    % computeMass computes the diagonal of the mass matrix
    %   computeMass( mesh, rho ) returns the lumped node diagonal of the mass 
    %   matrix for the given mesh

    m = zeros(mesh2d.N * 2, 1);
    el = mesh2d.el;
    for i = 1:size(el, 1)
        rho = mesh2d.materials(mesh2d.materialIndex(i)).rho;
        for j = 1:3
            %mass is the density multiplied by the area divided by the
            %number of nodes of the element
            m(el(i).t(j) * 2 - 1) = m(el(i).t(j) * 2 - 1) + (1/3) * el(i).area * rho;
            m(el(i).t(j) * 2)     = m(el(i).t(j) * 2)     + (1/3) * el(i).area * rho;
        end
    end
end