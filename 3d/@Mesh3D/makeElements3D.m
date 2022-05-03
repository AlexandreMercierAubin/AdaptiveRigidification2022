function [el] = makeElements3D(V, T)
% MAKEELEMENTS Creates elements from tetrahedrons.
%   makeElements( V, T ) returns an array of structures given 3D
%   points V (stored as rows), and tetrahedrons T defined by indices
%   (1 indexed and stored as rows).

el = repmat(struct('T', 1), size(T, 1), 1);
if size(T,2) == 3
    area = triangle_area(V,T);
end

for index = 1:size(T, 1)
    i = T(index, 1);
    j = T(index, 2);
    k = T(index, 3);
    if size(T,2) == 4
        l = T(index, 4);
        el(index).restLength = [V(j, :)' - V(i, :)', V(k, :)' - V(i, :)', V(l, :)' - V(i, :)'];
    else
        x = V(j, :)' - V(i, :)';
        y = V(k, :)' - V(i, :)';
        el(index).restLength = [x, y];
        n = cross(x,y);
        el(index).referenceCrossLength = norm(n);
        el(index).referenceSpaceNormal = n./norm(n);
    end
    
    el(index).t = T(index, :);
    
    if size(T,2) == 4
        area = abs(det(el(index).restLength) / 6);
        el(index).area = area;
    else
        el(index).area = area(i);
    end
end