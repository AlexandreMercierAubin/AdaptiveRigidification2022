function mesh = eleNodeLoader(filename, attributes, materials, scale)
    [V,I] = readNODE('3d/data/'+filename+'.node');
    [T,A] = readELE('3d/data/'+filename+'.ele');
    [F,J,K] = boundary_faces(T);
%     tetramesh(T,V);
%     vol = volume(V,T);
    if nargin < 4
        scale = [1,1,1];
    end
    V(:,1) = scale(1)*V(:,1);
    V(:,2) = scale(2)*V(:,2);
    V(:,3) = scale(3)*V(:,3);
    mesh = Mesh3D(V,T, attributes, materials, F, J);
end