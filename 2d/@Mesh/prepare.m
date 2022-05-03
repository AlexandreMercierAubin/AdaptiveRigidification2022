function prepare(obj, infMassForPinned)
    %forces some pinned dofs to have an infinite mass
    if infMassForPinned == 1
        for index = obj.pinnedInds
            obj.mass(index * 2 - 1:index * 2) = Inf; 
            obj.M(index * 2 - 1, index * 2 - 1) = Inf;
            obj.M(index * 2, index * 2) = Inf;
        end
    end
end