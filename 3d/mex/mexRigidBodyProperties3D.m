% mexRigidBodyProperties3D - Compute properties of 3D rigid bodies
%
% Computes the properties of the rigid bodies from scratch.  Note that 
% this is linear and probably pretty fast, but we could potentailly make
% this even faster with incremental updates to the rigid bodies (e.g., if 
% just a few elements and nodes join an existing rigid, the update is very
% small and fast).  Unfortunately, bookkeeping for such things seems a bit
% tricky.
% 
% [ com, comdot, mass, J, omega ] = mexRigidBodyProperties3D( numRigid, rigidIDbyVert, p, pdot, mass ) 
%
%  inputs:
%
%      numRigid        scalar, number of rigids (possibly zero)
%      rigidIDbyVert   #Vert by 1 vector of rigid IDs
%      p               #vertx 3 by 1 vector of positions
%      pdot            #vertx 3 by 1 vector of velocities
%      mass            #vertx 3 by 1 lumped diagonal mass matrix
%
%  outputs:
%
%      com                 3 by #rigid position of center of mass of rigid
%      comdot              3 by #rigid velocity of center of mass of rigid
%      mass                1 by #rigid mass of a rigid
%      rotMass             9 by #rigid rotational inertia of rigid
%      angularMomentum     3 by #rigid current angular momentum, don't forget to solve AngularVelocity = Inertia \ angularMomentum
%      VertexDisp          position of vertex wrt its COM
%      
% To compile type: mex -R2018a mexRigidBodyProperties3D.cpp