/*==========================================================
 * mexPGS2DwithJAinvJT.cpp - mex projected gauss seidel
 *
 * Solves friction problem J A^{-1} J' lambda = b problem with given 
 * friction bounds. Output is lambda.
 *
 * The calling syntax is:
 *
 * [ lambda ] = mexPGS2DwithJAinvJT( iterations, lambda, JAinvJT, b, mu, compliance );
 *
 * To compile type: mex -R2018a mexPGS2DwithJAinvJT.cpp
 *========================================================*/

#include "mex.h"
#include "blas.h"
#include <algorithm>
#include <cmath>
/* 
 * The gateway function
 * ( iterations, lambda, JAinvJT, b, mu, compliance )
 * For N contacts, most vectors will have 2N entires
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    double *lambda;             /* 2nContacts x 1 (output, copied from input) */
    double *b;                  /* 2nContacts x 1 */
    double *mu;                 /* 1 x nContacts */
    double compliance;          /* scalar, but could be per contact?? */
    double *JAinvJT;            /* 2nContacts x 2nContacts */    
    size_t nContacts;           /* number of contacts */
    
    /* check for proper number of arguments */
    if ( nrhs != 6 ) {
        mexErrMsgIdAndTxt("ARP:mexPGS2DwithJAinvJT:nrhs","Six inputs required.");
    }
    if ( nlhs != 1 ) {
        mexErrMsgIdAndTxt("ARP:mexPGS2DwithJAinvJT:nlhs","One outputs required.");
    }
    /* make sure the first input argument is scalar */
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("ARP:mexPGS2DwithJAinvJT:notScalar","iterations must be a scalar.");
    }
    /* make sure the second input argument is type double */
    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ) {
        mexErrMsgIdAndTxt("ARP:mexPGS2DwithJAinvJT:notDouble","lambda vector must be type double.");
    }
    /* make sure the compliance input argument is scalar */
    if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || mxGetNumberOfElements(prhs[5])!=1 ) {
        mexErrMsgIdAndTxt("ARP:mexPGS2DwithJAinvJT:notScalar","compliance must be a scalar.");
    }

    /* TODO: more input checking!!! or play dangerous and assume perfect input? */
 
    int iterations = (int) mxGetScalar( prhs[0] );
    compliance = mxGetScalar( prhs[5] );
    b = mxGetDoubles( prhs[3] );
    mu = mxGetDoubles( prhs[4] );
    nContacts = mxGetN( prhs[4] );              // get number of contacts, and safer not to trust sizes
    size_t nContactsx2 = nContacts*2;           // get size of lambda, and safer not to trust sizes
    plhs[0] = mxDuplicateArray( prhs[1] );      // create the output vector by making a deep copy
    lambda = mxGetDoubles( plhs[0] );           // get the output vector data so we can modify it
    JAinvJT = mxGetDoubles( prhs[2] );          // get dense matrix values
    
    double Jcdv;
    int i;               // index for lambda
    double Dii;          // diagonal of the matrix

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
            i = contact_index*2;
            
            Jcdv = 0;
            for ( size_t col = 0; col < nContactsx2; col++ )  {
                Jcdv += JAinvJT[ i + col*nContactsx2 ] * lambda[ col ]; 
            }
            Dii = JAinvJT[ i + i*nContactsx2 ];
            lambda[i] = ( lambda[i]*Dii - b[i] - Jcdv ) / ( Dii + compliance );
            if (std::isnan(lambda[i])) {
                lambda[i] = 0;
            }
            lambda[i] = std::max(lambda[i], 0.0);

            // second, the tangent
            // if no friction do nothing
            if ( mu[contact_index] == 0 ) {
                continue;
            }
            i++;
            double ub = mu[contact_index] * lambda[i-1];
            Jcdv = 0;
            for ( size_t col = 0; col < nContactsx2; col++ )  {
                Jcdv += JAinvJT[ i + col*nContactsx2 ] * lambda[ col ]; 
            }
            Dii = JAinvJT[ i + i*nContactsx2 ];
            lambda[i] = ( lambda[i]*Dii - b[i] - Jcdv ) / Dii; // no compliance for friction!
            if (std::isnan(lambda[i])) {
                lambda[i] = 0;
            }
            lambda[i] = std::min( std::max(lambda[i], -ub), ub);
        }
    }
}
