% mexPGS2DwithJAinvJT.cpp - mex projected gauss seidel
%
% Solves friction problem J A^{-1} J' lambda = b problem with given 
% friction bounds. Output is lambda.
%
% The calling syntax is:
%
% [ lambda ] = mexPGS2DwithJAinvJT( iterations, lambda, JAinvJT, b, mu, compliance );
%
%   iterations         scalar
%   lambda             2nContacts x 1 
%   JAinvJT            2nContacts x 2nContacts dense matrix
%   b                  2nContacts x 1 
%   mu                 1 x nContacts,  we might not care if row or col
%   compliance         scalar, but could later become one per contact 
%
% To compile type: mex -R2018a mexPGS2DwithJAinvJT.cpp
