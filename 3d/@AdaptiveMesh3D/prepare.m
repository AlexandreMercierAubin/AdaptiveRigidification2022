function prepare(obj, infMassForPinned)
    prepare@Mesh(obj, infMassForPinned);
    obj.updateRigidState();
end