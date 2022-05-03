% mexSTVK3D.cpp - computes negative area weighted STVK dpsidF and d2psidF2
%
% [ ii, jj, Cvals, dpsidF ] = mexSTVK3D( F, V, mu, lambda );
%
% Input:
%  F        9x#E deformation gradient of each element (or a vector of size 9x#E)
%  V        #E volume of each element
%  mu       #E Lamé parameter of each element
%  lambda   #E Lamé parameter of each element
% 
% Output:
%  ii         81x#E by 1 row indices for each C block in a sparse matrix
%  jj         81x#E by 1 col indices for each C block in a sparse matrix
%  Cvals      81x#E by 1 energy Hessian wrt F for each eleemnt (scaled by negative volume)
%  dpsidF     9x#E by 1 energy gradient (scaled by negative volume)
%
% To compile type: mex -R2018a mexSTVK3D.cpp
%
% Notes:
%
% Code generated by scratch/codeGenSTVKHessian3D.m
%
% It is perhaps a tiny bit slower to build the sparse matrix with the 
% matlab call sparse( ii, jj, Cvals ), rather than creating it here, but 
% it could also be more useful in the long run to have the non-sparse 
% version of the Cvals (e.g., a custom non-assembled matrix multiply for
% the quick solve conjugate gradient).