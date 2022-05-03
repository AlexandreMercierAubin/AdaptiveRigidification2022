function B = computeB3D(V,T)
    % B = computeB(V,T)
    % computes the kinematic mapping B.
    % V #vertices X 3 entries of x,y,z coordinates of vertices.
    % T #elements X 4 entries of vertex numbers tied to elements.
    B = sparse(9*size(T,1), 3*max(T(:)));
    for i = 1:size(T,1)
        %Fetching the position of eatch vertex
        t = T(i,:);
        p0 = V(t(1),:);
        p1 = V(t(2),:);
        p2 = V(t(3),:);
        p3 = V(t(4),:);
        
        %femdefo's Dm matrix, this is essentially a material frame of the
        %undeformed mesh
        materialFrame = [p1-p0;p2-p0;p3-p0]';
        materialFrameInv = inv(materialFrame);
        
        %local B block for individual elements
        D = [-ones(1,3)*materialFrameInv;materialFrameInv] ;
        localB = zeros(9,12);
        localB(1:3,1:3:end) = D';
        localB(4:6,2:3:end) = D';
        localB(7:9,3:3:end) = D';
        
        %plug the local B into a global B matrix that cont
        rowEnd = i*9;
        cols = reshape([t*3-2;t*3-1;t*3],1,[]);
        B(rowEnd - 8: rowEnd,cols) = localB;
    end

end

