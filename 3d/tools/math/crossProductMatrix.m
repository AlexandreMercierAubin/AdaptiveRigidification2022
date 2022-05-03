function [crossMat] = crossProductMatrix(m)
% makes a cross product matrix from a vector omega
crossMat = [0,   -m(3),   m(2);
            m(3),    0,  -m(1);
           -m(2) ,m(1),     0];
end

