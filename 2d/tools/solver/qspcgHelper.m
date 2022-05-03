function Av = qspcgHelper( h, M, Md, bigB, bigAlpha1, bigC, v )
% QSPCGHELPER Provides a multiplication by system matrix A without assembly
% and does so with an elimination of common subexpressions (something that
% cannot be done with an anonymous function).
    CBv = bigC*(bigB*v);
    tmp = h*(bigAlpha1*CBv) + h^2*CBv;
    Av = M*v + h*(Md*v) - bigB'*tmp;
end