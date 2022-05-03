function [] = genHandMadePOLY()
% GENHANDMADEPOLY Creates small examples programmatically and saves to POLY
% ELE and NODE FILES.
 
%     [p,t] = makeSingleTri();
%     writePT2POLY( 'data/singleTri', p, t );
 
%     [p,t] = makeTwoTri();
%     writePT2POLY( 'data/twoTri', p, t );
     
%    [p,t] = makeTwoTriSymmetric();
%    writePT2POLY( 'data/twoTriSym', p, t );
    
     [p,t] = makeFourTri();
     writePT2POLY( 'data/FourTri', p, t );

     [p,t] = makeFourTri();
     writePT2POLY( 'data/SixTri', p, t );

    
    function [p,t] = makeSingleTri()
        p(1, 1) = 0;
        p(1, 2) = 0;
        p(2, 1) = 1;
        p(2, 2) = 0;
        p(3, 1) = 0;
        p(3, 2) = 1;
        t(1, 1) = 1;
        t(1, 2) = 2;
        t(1, 3) = 3;
    end

    function [p,t] = makeTwoTri()
        p(1, 1) = -1;
        p(1, 2) = -1;
        p(2, 1) = 1;
        p(2, 2) = -1;
        p(3, 1) = 1;
        p(3, 2) = 1;
        p(4, 1) = 2;
        p(4, 2) = -1;
        t(1, 1) = 1;
        t(1, 2) = 2;
        t(1, 3) = 3;
        t(2, 1) = 3;
        t(2, 2) = 2;
        t(2, 3) = 4;
    end

    function [p,t] = makeTwoTriSymmetric()
        p = [ -2  0;
               0 -1;
               0  1;
               2  0 ];
        t = [ 1 2 3;
              3 2 4 ];
    end

    function [p,t] = makeFourTri()
        % Strangely, when the two rigids on the end have the same mass the
        % momentum is conserved.  If they have different masses then
        % momentum is not conserved.
        p = [ -1  0;
               0 -1;
               0  1;
               1  0;
              -1 -1;
               1.25  1.25; %0.75;
               0  0 ];
        t = [ 1 2 7;
              2 4 7;
              4 3 7;
              3 1 7;
              5 2 1;
              6 3 4 ];
    end

    function [p,t] = makeSixTri()
        p = [ -1  0;
               0 -1;
               0  1;
               1  0;
              -1 -1;
               1  0.5 ];
        t = [ 1 2 3;
              3 2 4;
              5 2 1;
              6 3 4 ];
    end



    function writePT2POLY( filenameroot, P, T )
        mesh = Mesh( P, T, 1, 1, 1, 0, 0 ); % create a mesh just for the boundary edges    
        writeNODE( sprintf('%s.node', filenameroot), P );
        writeELE( sprintf('%s.ele', filenameroot), T );  % doesn't support attributes :(
        % POLY file isn't really even used.... but hey!
        edges = vertcat( mesh.boundaryEdges{:} );
        writePOLY_triangle( sprintf('%s.poly', filenameroot), P, edges, [] ); % no holes, and not currently needed, but pehraps later?
    end

end

