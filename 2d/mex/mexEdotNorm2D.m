% mexEdotNorm2D.cpp - computes 2D strain rate Frobeneneous norm squared
%
% [ eDotNormsSquared ] = mexEdotNorm2D( F, Fdot );
%
% Input:
%  F           4x#E deformation gradient of each element (or a vector of size 4x#E)
%  Fdot        4x#E deformation gradient velocity of each element (or a vector of size 4x#E)
%
% Output:
%  EdotNormSquared        #E Frobeneus norm squared of strain rate
%
% To compile type: mex -R2018a mexEdotNorm2D.cpp