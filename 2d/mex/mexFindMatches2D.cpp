/*==========================================================
 * mexFindMatches2D.cpp - mex find matches for warmstarting
 *
 * Find the matches between new contacts and old contacts. 
 *
 * The calling syntax is:
 *
 * [ warmStartLambdas, isNewContact ] = mexFindMatches2D( oldIDs, oldLambdas, newIDs );
 *
 * oldIDs       number of contacts by 1 array (possibly empty) of long
 * oldLambdas   2xcontacts by 1 array (possibly empty) of double
 * newIDs       number of new contacts y 1 array, non-empty
 * 
 * NOTE:
 *
 * If nothing changes between frames, it might be worthwhile to check if
 * all the contacts match, and if they do, simply copy the lambda.  Or copy 
 * those we can at the beginning, and use the map for the rest?
 *
 * To compile type: mex -R2018a mFindMatches.cpp
 *========================================================*/

#include "mex.h"
#include "blas.h"

#include <unordered_map> 

/* 
 * The gateway function
 * ( oldIDs, oldLambdas, newIDs )
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

     /* check for proper number of arguments */
    if ( nrhs != 3 ) {
        mexErrMsgIdAndTxt("ARP:mFindMatches:nrhs","3 inputs required.");
    }
    if ( nlhs != 2 ) {
        mexErrMsgIdAndTxt("ARP:mFindMatches:nlhs","Two outputs required.");
    }
    
    int64_T *oldIDs = (int64_T *)mxGetData( prhs[0] );
    size_t nOld = mxGetNumberOfElements( prhs[0] );

    double *oldLambdas = mxGetDoubles( prhs[1] );
    
    int64_T *newIDs = (int64_T *)mxGetData( prhs[2] );
    size_t nNew = mxGetNumberOfElements( prhs[2] );
    
    plhs[0] = mxCreateDoubleMatrix( nNew*2, 1, mxREAL ); // will be initialized to zero
    double *warmStartLambdas = mxGetDoubles( plhs[0] );
    plhs[1] = mxCreateLogicalMatrix( nNew, 1 );
    mxLogical *isNewContact = mxGetLogicals( plhs[1] );
    
    std::unordered_map<int64_T,int> umap;
    umap.reserve( nOld*2 );
    
    for ( int i = 0; i < nOld; i++ ) {
        umap[ oldIDs[i] ] = i;
    }
    for ( int i = 0; i < nNew; i++ ) {
        auto entry = umap.find( newIDs[i] );
        if ( entry == umap.end() ) {
            isNewContact[i] = true;
        } else {
            isNewContact[i] = false;
            int j = entry->second;
            warmStartLambdas[ i*2 ] = oldLambdas[ j*2 ];
            warmStartLambdas[ i*2+1 ] = oldLambdas[ j*2+1 ];
        }
    }
}