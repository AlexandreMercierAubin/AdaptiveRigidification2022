% mexFindMatches2D.cpp - mex find matches for warmstarting
%
% Find the matches between new contacts and old contacts. 
%
% The calling syntax is:
%
% [ warmStartLambdas, isNewContact ] = mexFindMatches2D( oldIDs, oldLambdas, newIDs );
%
% oldIDs       number of contacts by 1 array (possibly empty) of long
% oldLambdas   2xcontacts by 1 array (possibly empty) of double
% newIDs       number of new contacts y 1 array, non-empty
% 
% NOTE:
%
% If nothing changes between frames, it might be worthwhile to check if
% all the contacts match, and if they do, simply copy the lambda.  Or copy 
% those we can at the beginning, and use the map for the rest?
%
% To compile type: mex -R2018a mFindMatches.cpp
