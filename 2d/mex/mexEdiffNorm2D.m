% mexEdiffNorm2D.cpp - computes 2D strain difference Frobeneneous norm squared
%
% [ eDotNormsSquared ] = mexEdotNorm2D( Fa, Fb );
%
% Input:
%  Fa           4x#E deformation gradient of each element (or a vector of size 4x#E)
%  Fb           4x#E deformation gradient of each element (or a vector of size 4x#E)
%  h            timestep
%
% Output:
%  EdiffNormSquared        #E Frobeneus norm squared of (strain difference
%                             divided by h), i.e., an approximation of strain rate
%
% To compile type: mex -R2018a mexEdiffNorm2D.cpp