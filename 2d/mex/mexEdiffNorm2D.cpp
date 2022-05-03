/*==========================================================
 * mexEdiffNorm.cpp - computes 2D strain difference Frobeneneous norm squared
 *
 * [ eDotNormsSquared ] = mexEdotNorm( Fa, Fb, h );
 *
 * Input:
 *  Fa      4x#E deformation gradient of each element (or a vector of size 4x#E) at one time step
 *  Fb      4x#E deformation gradient of each element (or a vector of size 4x#E) at another time step
 *  h       timestep
 *
 * Output:
 *  EdiffNormSquared     #E Frobeneus norm squared of strain difference 
 *                          divided by h), i.e., an approximation of strain rate
 *
 * To compile type: mex -R2018a mexEdiffNorm.cpp
 *========================================================*/

#include "mex.h"
#include "blas.h"
#include <math.h>

/* 
 * The gateway function
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    double *Fa = mxGetDoubles(prhs[0]);
    size_t sizeF = mxGetM(prhs[0]);
    double *Fb = mxGetDoubles(prhs[1]);
    double h = mxGetScalar( prhs[2] );

    int sizeEdiffNorm = sizeF/4.0;
    plhs[0] = mxCreateDoubleMatrix( sizeEdiffNorm, 1, mxREAL );
    double *EdiffNormSquared = mxGetDoubles(plhs[0]);
    
    for(size_t i = 0; i < sizeEdiffNorm; ++i) {
        int pos = i*4;
        double Fa1_1 = Fa[pos];
        double Fa2_1 = Fa[pos+1];
        double Fa1_2 = Fa[pos+2];
        double Fa2_2 = Fa[pos+3];
        
        double Fb1_1 = Fb[pos];
        double Fb2_1 = Fb[pos+1];
        double Fb1_2 = Fb[pos+2];
        double Fb2_2 = Fb[pos+3]; 
        
        double a = Fa1_1*Fa1_1 + Fa2_1*Fa2_1 - Fb1_1*Fb1_1 - Fb2_1*Fb2_1;
        double b = Fa1_1*Fa1_2 + Fa2_1*Fa2_2 - Fb1_1*Fb1_2 - Fb2_1*Fb2_2;
        double c = Fa1_2*Fa1_2 + Fa2_2*Fa2_2 - Fb1_2*Fb1_2 - Fb2_2*Fb2_2;

        EdiffNormSquared[i] = (0.25*a*a + 0.5*b*b + 0.25*c*c)/h/h;
    }
}
