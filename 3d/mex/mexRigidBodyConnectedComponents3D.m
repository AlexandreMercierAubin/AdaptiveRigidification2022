% mexRigidBodyConnectedComponents3D - Find rigid components for 3D mesh
%
% Computes connected rigid components based on Edot history, but likewise 
% prevents rigidification in any elements that are currently showing a 
% large approxEdot from the quick solve.  This method uses the triangle
% connectivity graph (possibly just edge adjacency, but could also be one
% of two other alternatives to avoid hinges).
%
% The general strategy does a BFS starting form the towPinTris elements, 
% but then simply walking over all tris afterward and restarting the BFS
% for a rigid component at any elements that have not yet been visited.
% Temp memory is needed to keep track of vertices and triangles visited.
% Rigidificaiton rules will combine the input EdotHistory, EdotApprox, 
% along with one of the following:
% 1) Tri-adj across edges, and a first come first serve greedy approach to
% giving a vertex to a rigid component. 
% 2) Tri-adj across edges and vertices, recognizing that hinges will be 
% welded
% 3) a more complex multi pass approach of identiyfing vertices that have 
% all triangles ready to rigidfy, then allowing triangles that have all 
% verts ready to go become rigid (this way tri-adj across edges and 
% vertices will weld, but all surrounding triangles will be in a low Edot
% state so it will be 100% Kosher!
% 
% [ numRigid, rigidIDbyTri rigidIDbyVert, isElasticVert, isBoundaryVert, isElasticTri ] = mexRigidConnectedComponents3D( ... 
%     EdotHistory, EdotApprox, adjMatrix, triVertIDint32, numVerts, pinnedVerts ) 
%
%  inputs:
%
%      EdotHistory     #Tri by 3? Edot history of elements
%      EdotApprox      #tri by 1 result from quick solve no constraints
%      adjMatrix       #tri by #tri SPARSE symmetric matrix of triangle adj
%      triVertIDint32  #tri by 3 (4 for tets) indices (starting from 1) of verts of each element
%      numVerts        scalar number of vertices
%      ispinnedVerts     #vert by 1 logical identifying pinned vertices
%      verticeValance  #vert by 1 array of the valances
%
%  outputs:
%
%      numRigid        scalar, number of rigids (possibly zero)
%      rigidIDbyTri	#Tri by 1 vector of rigid IDs 
%      rigidIDbyVert   #Vert by 1 vector of rigid IDs
%      isVertElastic   #Vert by 1 logical identifying rigid verts (i.e., rigidIDbyVert == -1, and avoid a find on the matlab side)
%      isVertBoundary  #Vert by 1 logical identifying rigid verts on the boundary of a rigid
%      isTriElastic    #Tri by 1 logical identifying rigid tris (i.e., rigidIDbyTri == -1, and avoid a find on the matlab side) 
%
%
% To compile type: mex -R2018a mexRigidBodyConnectedComponents3D.cpp
