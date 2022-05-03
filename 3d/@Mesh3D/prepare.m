function prepare(obj, infMassForPinned)
    if infMassForPinned == 1
        for index = obj.pinnedInds
            obj.mass(index * 3 - 2:index * 3) = Inf; 
            obj.M(index * 3 - 2, index * 3 - 2) = Inf;
            obj.M(index * 3 - 1, index * 3 - 1) = Inf;
            obj.M(index * 3, index * 3) = Inf;
        end
    end
end