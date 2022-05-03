/*==========================================================
 * mexEdotNorm.cpp - computes 2D strain rate Frobeneneous norm squared
 *
 * [ eDotNormsSquared ] = mexEdotNorm( F, Fdot );
 *
 * Input:
 *  F           4x#E deformation gradient of each element (or a vector of size 4x#E)
 *  Fdot        4x#E deformation gradient velocity of each element (or a vector of size 4x#E)
 *
 * Output:
 *  EdotNormSquared     #E Frobeneus norm squared of strain rate
 *
 * To compile type: mex -R2018a mexEdotNorm.cpp
 *========================================================*/

#include "mex.h"
#include "blas.h"
#include <math.h>

/* 
 * The gateway function
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    double *F = mxGetDoubles(prhs[0]);
    size_t sizeF = mxGetM(prhs[0]);
    double *FDot = mxGetDoubles(prhs[1]);
    int sizeEdotNorm = sizeF/4.0;
    plhs[0] = mxCreateDoubleMatrix( sizeEdotNorm, 1, mxREAL );
    double *EdotNormSquared = mxGetDoubles(plhs[0]);
    
    for(size_t i = 0; i < sizeEdotNorm; ++i){
        int pos = i*4;
        double F1_1 = F[pos];
        double F2_1 = F[pos+1];
        double F1_2 = F[pos+2];
        double F2_2 = F[pos+3];
        
        double Fd1_1 = FDot[pos];
        double Fd2_1 = FDot[pos+1];
        double Fd1_2 = FDot[pos+2];
        double Fd2_2 = FDot[pos+3]; 
        
        double a = F1_1*Fd1_2 + F1_2*Fd1_1 + F2_1*Fd2_2 + F2_2*Fd2_1;
        double b = F1_1*Fd1_1 + F2_1*Fd2_1;
        double c = F1_2*Fd1_2 + F2_2*Fd2_2;

        EdotNormSquared[i] = 0.5*a*a + b*b + c*c;
        //EdotNormSquared[i] = 0.25*a*a + 0.5*b*b + 0.25*c*c;
    }
}
