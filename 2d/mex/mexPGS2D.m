% mexPGS2D.c - mex projected gauss seidel
%
% Solves friction problem J A^{-1} J' lambda = b problem with given 
% friction bounds.  Note T = A^{-1} J', and Dii contains the diagonals of
% the left hand side matrix.  Outputs are lambda and deltav, while 
% warmstart values of lambda and deltav are passed as arguments. 
%
% The calling syntax is:
%
% [ lambda, deltav ] = mexPGS2D( iterations, lambda, deltav, T, Dii, b, JcT, mu, compliance );
%
%   iterations         scalar
%   lambda             2nContacts x 1 
%   deltav             nDOFs x 1
%   T                  nDOFs x 2nContacts
%   Dii                1 x 2nContacts, we might not care if row or col
%   b                  2nContacts x 1 
%   JcT                nDOFs x 2nContacts sparse matrix
%   mu                 1 x nContacts,  we might not care if row or col
%   compliance         scalar, but could later become one per contact 
%
% To compile type: mex -R2018a mexPGS2D.cpp
