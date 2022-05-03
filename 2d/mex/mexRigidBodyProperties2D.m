% mexRigidBodyProperties2D - mex rigid body property computation
%
% Computes the properties of the rigid bodies from scratch.  Note that 
% this is linear and probably pretty fast, but we could potentailly make
% this even faster with incremental updates to the rigid bodies (e.g., if 
% just a few elements and nodes join an existing rigid, the update is very
% small and fast).  Unfortunately, bookkeeping for such things seems a bit
% tricky.
% 
% [ com, comdot, mass, J, omega ] = mexRigidBodyProperties2D( numRigid, rigidIDbyVert, p, pdot, mass ) 
%
%  inputs:
%
%      numRigid         scalar, number of rigids (possibly zero)
%      rigidIDbyVert    #Vert by 1 vector of rigid IDs
%      p                #vertx2 by 1 vector of positions
%      pdot             #vertx2 by 1 vector of velocities
%      mass             #vertx2 by 1 lumped diagonal mass matrix
%
%  outputs:
%
%      com         2 by #rigid position of center of mass of rigid
%      comdot      2 by #rigid velocity of center of mass of rigid
%      mass        1 by #rigid mass of a rigid
%      rotMass     1 by #rigid rotational inertia of rigid (will be 3x3 matrix per rigid in 3D)
%      omega       1 by #rigid current angular velocity
%      VertexDisp  position of vertex wrt its COM
%      
% To compile type: mex -R2018a mexRigidBodyProperties2D.cpp