function prepare(obj, infMassForPinned)
    %forces some pinned dofs to have an infinite mass
    %initializes some rigid state in case we start rigidified
    prepare@Mesh(obj, infMassForPinned);
    obj.updateRigidState();
end