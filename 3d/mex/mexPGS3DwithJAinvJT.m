% mexPGS3DwithJAinvJT.cpp - mex projected gauss seidel
%
% Solves friction problem J A^{-1} J' lambda = b problem with given
% friction bounds.  Output is lambda, while warmstart values of lambda are
% passed as arguments.
%
% The calling syntax is:
%
% [ lambda ] = mexPGS3DwithJAinvJT( iterations, lambda, JAinvJT, b, mu, compliance );
%
%   iterations         scalar
%   lambda             3nContacts x 1 
%   JAinvJT            3nContacts x 3nContacts dense matrix
%   b                  3nContacts x 1 
%   mu                 1 x nContacts,  we might not care if row or col
%   compliance         scalar, but could later become one per contact  
% 
% To compile type: mex -R2018a mexPGS3DwithJAinvJT.cpp