/*==========================================================
 * mexPGS3D.cpp - mex projected gauss seidel
 *
 * Solves friction problem J A^{-1} J' lambda = b problem with given 
 * friction bounds.  Note T = A^{-1} J', and Dii are the diagonals of
 * the left hand side matrix.  Outputs are lambda and deltav, while 
 * warmstart values of lambda and deltav are passed as arguments. 
 *
 * The calling syntax is:
 *
 * [ lambda, deltav ] = mexPGS3D( iterations, lambda, deltav, T, Dii, b, JcT, mu, compliance );
 *
 * Notice the Jacobian transpose must be provided because of the compressed column format!  
 *
 * To compile type: mex -R2018a mexPGS3D.cpp
 *========================================================*/

#include "mex.h"
#include "blas.h"
#include <algorithm>
#include <cmath>
/* 
 * The gateway function
 * ( iterations, lambda, deltav, T, Dii, b, JcT, mu, compliance )
 * JcT must be sparse
 * For N contacts, many vectors will have 3N entires
 * Size of T is M-by-3N
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    double *lambda;             /* 3nContacts x 1 (output, duplicated from input) */
    double *deltav;             /* nDOFs x 1  (output, duplicated from input) */
    double *b;                  /* 3nContacts x 1 */
    double *Dii;                /* 1 x 3nContacts, do we care row/col? */
    double *mu;                 /* 1 x nContacts */
    double compliance;          /* scalar, but could be per contact?? */
    double *T;                  /* nDOFs x 3nContacts */
    
    size_t nContacts;           /* number of contacts */
    size_t nDOFs;               /* number of dofs */
    
    /* check for proper number of arguments */
    if ( nrhs != 9 ) {
        mexErrMsgIdAndTxt("ARP:mpgs:nrhs","Nine inputs required.");
    }
    if ( nlhs != 2 ) {
        mexErrMsgIdAndTxt("ARP:mpgs:nlhs","Two outputs required.");
    }
    /* make sure the first input argument is scalar */
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("ARP:mpgs:notScalar","iterations must be a scalar.");
    }
    /* make sure the second input argument is type double */
    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ) {
        mexErrMsgIdAndTxt("ARP:mpgs:notDouble","lambda vector must be type double.");
    }

    /* make sure the first input argument is scalar */
    if( !mxIsDouble(prhs[8]) || mxIsComplex(prhs[8]) || mxGetNumberOfElements(prhs[8])!=1 ) {
        mexErrMsgIdAndTxt("ARP:mpgs:notScalar","compliance must be a scalar.");
    }

    /* TODO: more input checking!!! or play dangerous and assume perfect input? or is this stuff slow?*/
 
    /* get the value of the scalar input  */
    int iterations = (int) mxGetScalar( prhs[0] );
    compliance = mxGetScalar( prhs[8] );

    /* create a pointer to the data in the input matrix  */
    T = mxGetDoubles( prhs[3] );
    Dii = mxGetDoubles( prhs[4] );
    b = mxGetDoubles( prhs[5] );
    mu = mxGetDoubles( prhs[7] );

    /* get dimensions of the input matrix */
    nContacts = mxGetN( prhs[7] ); // safer not to trust sizes
    nDOFs = mxGetM( prhs[2] );

    /* create the output vectors by making deep copies*/
    plhs[0] = mxDuplicateArray( prhs[1] );
    plhs[1] = mxDuplicateArray( prhs[2] );

    /* get the output vector data so we can modify it */
    lambda = mxGetDoubles( plhs[0] );
    deltav = mxGetDoubles( plhs[1] );
    
    
    /* call the computational routine */
//    arrayProduct(multiplier,inMatrix,outMatrix,(mwSize)ncols);
    
    //analyze_sparse( prhs[6] ) ;
    
    double *JcTvals = mxGetDoubles( prhs[6] );
    mwIndex *ir = mxGetIr( prhs[6] );
    mwIndex *jc = mxGetJc( prhs[6] );

    // for sparse matrix vector product, Jc(i,:)*deltav, see:
    //    https://itectec.com/matlab/matlab-a-sparse-matrix-vector-multiplication-in-c-mex-programming/
    // see also
    // https://www.mathworks.com/matlabcentral/answers/402581-how-does-matlab-internally-store-sparse-arrays-calling-mkl-s-sparse-blas-library-for-sparse-matrix
    
    
    
    // 99% sure that pointer order is row then column
    //    double *T;                  /* nDOFs x 2nContacts */
//     for ( int i = 0; i < nDOFs*2*nContacts; i++ ) {
//         mexPrintf( "%g\n", T[i] );
//     }
//     return;
    
    
    double oldLambda;
    double Jcdv;
    size_t starting_row_index;
    size_t stopping_row_index;
    double deltaLambda;         
    size_t offset;                 // offset to T column data
    int i;                      // index for lambda
    
    for ( int iter = 0; iter < iterations; iter++ ) {
        // mexPrintf( "Iteration = %d\n",  i );
        
        // This index is for walking through the JcT values, and note that
        // we will use up values in the first column for the normal, and
        // then use up values in the next column for the tangent, and then
        // continue to the next contact repeating.  Thus this index is set
        // to zero outside of the loop over all contacts.
        int index = 0;  
        for ( int contact_index = 0; contact_index < nContacts; contact_index++ ) {

            // first normal, then tangent
            i = contact_index*3;
            
            // Compute Jcdv as dot product of column i of JcT with dv
            starting_row_index = jc[i];
            stopping_row_index = jc[i+1];
            if ( starting_row_index == stopping_row_index ) {
                mexPrintf( "zero constraint row?  quitting" );
                return;
            }
            Jcdv = 0;
            for ( size_t current_row_index = starting_row_index; current_row_index < stopping_row_index; current_row_index++ )  {
                size_t row = ir[current_row_index];
                Jcdv += JcTvals[index++] * deltav[row];
            }
            
            // Compute lambda update
            oldLambda = lambda[i];
            lambda[i] = ( lambda[i]*Dii[i] - b[i] - Jcdv ) / ( Dii[i] + compliance );
            if (std::isnan(lambda[i])) {
                lambda[i] = 0;
            }
            lambda[i] = std::max(lambda[i], 0.0);

            // Compute deltav update
            // recall double *T has size nDOFs x 3nContacts 
            deltaLambda = lambda[i] - oldLambda;
            offset = i * nDOFs; // to index values in the desired column
            for ( int j = 0; j < nDOFs; j++ ) {
                deltav[j] += T[offset+j] * deltaLambda;
            }

            //lambda of normal
            double ub = mu[contact_index] * lambda[i];
            // second, the tangent
            for (int k = 1; k < 3; ++k){
                i++;

                // Might not be that common, but if lambda[i] is zero, and was
                // previously zero then we can skip the tangent, and otherwise
                // just need to go through once... that is, could also check 
                // that the normal force is zero but the tangent is not.
                // TODO: worry about the correctness of this!
                if ( mu[contact_index] == 0 ) { //|| lambda[i-1] == 0 && deltaLambda == 0 ) {
                    index += jc[i+1] - jc[i];
                    continue;
                }
                
                // Compute Jcdv as dot product of column i of JcT with dv
                Jcdv = 0;
                starting_row_index = jc[i];
                stopping_row_index = jc[i+1];
                for ( size_t current_row_index = starting_row_index; current_row_index < stopping_row_index; current_row_index++ )  {
                    size_t row = ir[current_row_index];
                    Jcdv += JcTvals[index++] * deltav[row];
                }

                // Compute lambda update
                oldLambda = lambda[i];
                lambda[i] = ( lambda[i]*Dii[i] - b[i] - Jcdv ) / Dii[i]; // no compliance for friction!
                if (std::isnan(lambda[i])) {
                    lambda[i] = 0;
                }
                lambda[i] = std::min( std::max(lambda[i], -ub), ub);

                // Compute deltav update
                // recall double *T has size nDOFs x 3nContacts 
                deltaLambda = lambda[i] - oldLambda;
                offset = i * nDOFs; // to index values in the desired column
                for ( int j = 0; j < nDOFs; j++ ) {
                    deltav[j] += T[offset+j] * deltaLambda;
                }
            }
        }
    }
}
