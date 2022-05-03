function [R U V] = polarDecomposition(F)
%POLDECOMP  Performs the polar decomposition of a regular square matrix.
%   Zoltán Csáti (2021). Polar decomposition (https://www.mathworks.com/matlabcentral/fileexchange/48735-polar-decomposition), MATLAB Central File Exchange. Retrieved March 17, 2021.
%   Alex: I added the normalization as suggested in the comments
%   [R U V] = POLDECOMP(F) factorizes a non-singular square matrix F such
%   that F=R*U and F=V*R, where
%   U and V are symmetric, positive definite matrices and
%   R is a rotational matrix
%
%   See also EIG, DIAG, REPMAT


% This kind of decomposition is often used in continuum mechanics so it is
% convenient to comment the code that way. From now, we use the matrix 
% formalism of tensors. C is the right Cauchy-Green deformation tensor, 
% F is the deformation tensor, lambda is the stretch.

% Check input
[m n] = size(F);
if m ~= n
    error('Matrix must be square.');
end

C = F'*F;
[Q0, lambdasquare] = eig(C);
Q0 = Q0/norm(Q0);
lambda = sqrt(diag((lambdasquare))); % extract the components
% Uinv is the inverse of U and is constructed with the help of Q0. Uinv is
% produced in the same base as F not in the base of its eigenvectors.
Uinv = repmat(1./lambda',size(F,1),1).*Q0*Q0';
% Using the definition, R, U and V can now be calculated
R = F*Uinv;
U = R'*F;
V = F*R';