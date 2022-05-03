/*==========================================================
 * mexBsplit.cpp
 *
 *  inputs:
 *
 *      Btransposed
 *      activeBrows
 *
 *  outputs:
 *
 *      ii1
 *      jj1
 *      C1
 *      ii2
 *      jj2
 *      C2
 *
 * To compile type: mex -R2018a mexBSplit.cpp
 *========================================================*/

#include "mex.h"
#include "blas.h"
#include <map>
#include <vector>
#include <queue>

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    double* Btransposed = mxGetDoubles(prhs[0]);   
    mwIndex* ir = mxGetIr( prhs[0] );
    mwIndex* jc = mxGetJc( prhs[0] );
    double *pr = mxGetPr (prhs[0]);
    size_t numcolsBt = mxGetN(prhs[0]);
    size_t numrowsBt = mxGetM(prhs[0]);
    
    double* activeBrows = mxGetDoubles(prhs[1]);
    size_t numActiveRow = mxGetM(prhs[1]);
    //TODO: add some asserts
    //TODO: change those for vectors and then use the size to init the double matrices
    
    // create the output vectors Belastic
    std::vector<double> ii1;
    std::vector<double> jj1;
    std::vector<double> C1; 

    // create the output vectors Brigid
    std::vector<double> ii2; 
    std::vector<double> jj2; 
    std::vector<double> C2;
    
    size_t counter = 0;
    size_t counterInactive = 0;
    for(size_t col = 0; col<numcolsBt; ++col){
        size_t starting_row_index = jc[col];
        size_t stopping_row_index = jc[col + 1];
        //remember, this is Btransposed
        bool isActiveRow = (counter < numActiveRow) && (activeBrows[counter] - 1 == col);
        if (isActiveRow) { counter++; }
        else { counterInactive++; }
        
        for ( size_t current_row_index = starting_row_index; current_row_index < stopping_row_index; current_row_index++ )  {
            size_t row = ir[current_row_index];
            double value = pr[current_row_index];
            //inverting rows and cols because the input is Btransposed
            if(isActiveRow){
                ii1.push_back((double)counter);
                jj1.push_back((double)row+1);
                C1.push_back((double)value);
            } else{
                ii2.push_back((double)counterInactive);
                jj2.push_back((double)row+1);
                C2.push_back((double)value);
            }
        }
    }
    
 	plhs[0] = mxCreateDoubleMatrix(C1.size(), 1, mxREAL );
    double* ii1mat = mxGetDoubles(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(C1.size(), 1, mxREAL );
    double* jj1mat = mxGetDoubles(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(C1.size(), 1, mxREAL );
    double* C1mat = mxGetDoubles(plhs[2]);
    for(size_t i = 0; i<C1.size();++i){
        ii1mat[i] = ii1[i];
        jj1mat[i] = jj1[i];
        C1mat[i] = C1[i];
    }
    
    // create the output vectors Brigid
    plhs[3] = mxCreateDoubleMatrix( C2.size(), 1, mxREAL );
    double* ii2mat = mxGetDoubles(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix( C2.size(), 1, mxREAL );
    double* jj2mat = mxGetDoubles(plhs[4]);
    plhs[5] = mxCreateDoubleMatrix( C2.size(), 1, mxREAL );
    double* C2mat = mxGetDoubles(plhs[5]);
    for(size_t i = 0; i<C2.size(); ++i){
        ii2mat[i] = ii2[i];
        jj2mat[i] = jj2[i];
        C2mat[i] = C2[i];
    }
}
