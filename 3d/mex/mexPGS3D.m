% mexPGS3D.cpp - mex projected gauss seidel
%
% Solves friction problem J A^{-1} J' lambda = b problem with given 
% friction bounds.  Note T = A^{-1} J', and Dii are the diagonals of
% the left hand side matrix.  Outputs are lambda and deltav, while 
% warmstart values of lambda and deltav are passed as arguments. 
%
% The calling syntax is:
%
% [ lambda, deltav ] = mexPGS3D( iterations, lambda, deltav, T, Dii, b, JcT, mu, compliance );
%
% Notice the Jacobian transpose must be provided because of the compressed column format!  
%
% To compile type: mex -R2018a mexPGS3D.cpp
