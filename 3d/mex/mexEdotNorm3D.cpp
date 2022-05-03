/*==========================================================
 * mexEdotNorm3D.cpp - computes 3D strain rate Frobeneneous norm squared
 *
 * [ eDotNormsSquared ] = mexEdotNorm3D( F, Fdot );
 *
 * Input:
 *  F           9x#E deformation gradient of each element (or a vector of size 4x#E)
 *  Fdot        9x#E deformation gradient velocity of each element (or a vector of size 4x#E)
 *
 * Output:
 *  EdotNormSquared        #E Frobeneus norm squared of strain rate
 *
 * To compile type: mex -R2018a mexEdotNorm3D.cpp
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
    double *Fd = mxGetDoubles(prhs[1]);
    int sizeEdot = sizeF/9;
    plhs[0] = mxCreateDoubleMatrix( sizeEdot, 1, mxREAL );
    double *EdotNormSquared = mxGetDoubles(plhs[0]);
    
    for(size_t i = 0; i < sizeEdot; ++i){
        int pos = i*9;
        double F1_1 = F[pos];
        double F2_1 = F[pos+1];
        double F3_1 = F[pos+2];
        double F1_2 = F[pos+3];
        double F2_2 = F[pos+4];
        double F3_2 = F[pos+5];
        double F1_3 = F[pos+6];
        double F2_3 = F[pos+7];
        double F3_3 = F[pos+8];

        double Fd1_1 = Fd[pos];
        double Fd2_1 = Fd[pos+1];
        double Fd3_1 = Fd[pos+2];
        double Fd1_2 = Fd[pos+3];
        double Fd2_2 = Fd[pos+4];
        double Fd3_2 = Fd[pos+5];
        double Fd1_3 = Fd[pos+6];
        double Fd2_3 = Fd[pos+7];
        double Fd3_3 = Fd[pos+8];

        double a = F1_1*Fd1_2 + F1_2*Fd1_1 + F2_1*Fd2_2 + F2_2*Fd2_1 + F3_1*Fd3_2 + F3_2*Fd3_1;
        double b = F1_1*Fd1_3 + F1_3*Fd1_1 + F2_1*Fd2_3 + F2_3*Fd2_1 + F3_1*Fd3_3 + F3_3*Fd3_1;
        double c = F1_2*Fd1_3 + F1_3*Fd1_2 + F2_2*Fd2_3 + F2_3*Fd2_2 + F3_2*Fd3_3 + F3_3*Fd3_2;
        double d = F1_1*Fd1_1 + F2_1*Fd2_1 + F3_1*Fd3_1;
        double e = F1_2*Fd1_2 + F2_2*Fd2_2 + F3_2*Fd3_2;
        double f = F1_3*Fd1_3 + F2_3*Fd2_3 + F3_3*Fd3_3;
        
        EdotNormSquared[i] = a*a/2 + b*b/2 + c*c/2 + d*d + e*e + f*f;
    }
}
