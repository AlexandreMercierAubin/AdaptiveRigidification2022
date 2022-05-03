% mexEdotNorm3D.cpp - computes 3D strain rate Frobeneneous norm squared
%
% [ eDotNormsSquared ] = mexEdotNorm3D( F, Fdot );
%
% Input:
%  F           9x#E deformation gradient of each element (or a vector of size 4x#E)
%  Fdot        9x#E deformation gradient velocity of each element (or a vector of size 4x#E)
%
% Output:
%  EdotNormSquared        #E Frobeneus norm squared of strain rate
%
% To compile type: mex -R2018a mexEdotNorm3D.cpp