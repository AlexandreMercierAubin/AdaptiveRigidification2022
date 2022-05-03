/*==========================================================
 * mRigidBodyProperties.c
 *
 * Computes the properties of the rigid bodies from scratch.  Note that 
 * this is linear and probably pretty fast, but we could potentailly make
 * this even faster with incremental updates to the rigid bodies (e.g., if 
 * just a few elements and nodes join an existing rigid, the update is very
 * small and fast).  Unfortunately, bookkeeping for such things seems a bit
 * tricky.
 * 
 * [ com, comdot, mass, J, omega ] = mRigidBodyProperties( numRigid, rigidIDbyVert, p, pdot, mass ) 
 *
 *  inputs:
 *
 *      numRigid        scalar, number of rigids (possibly zero)
 *      rigidIDbyVert   #Vert by 1 vector of rigid IDs
 *      p               #vertx2 by 1 vector of positions
 *      pdot            #vertx2 by 1 vector of velocities
 *      mass            #vertx2 by 1 lumped diagonal mass matrix
 *
 *  outputs:
 *
 *      com         2 by #rigid position of center of mass of rigid
 *      comdot      2 by #rigid velocity of center of mass of rigid
 *      mass        1 by #rigid mass of a rigid
 *      rotMass     1 by #rigid rotational inertia of rigid (will be 3x3 matrix per rigid in 3D)
 *      omega       1 by #rigid current angular velocity
 *      VertexDisp  position of vertex wrt its COM
 *      
 * To compile type: mex -R2018a mRigidBodyProperties.cpp
 *========================================================*/

#include "mex.h"
#include "blas.h"
#include <math.h>

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
        
    /* check for proper number of arguments */
    if ( nrhs != 5 ) {
        mexErrMsgIdAndTxt("ARP:mRigidBodyProperties:nrhs","5 inputs required.");
    }
    if ( nlhs != 6 ) {
        mexErrMsgIdAndTxt("ARP:mRigidBodyProperties:nlhs","6 outputs required.");
    }
 
    size_t numVert = mxGetNumberOfElements( prhs[1] );

    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("ARP:mRigidBodyProperties:size","numRigid must be a scalar.");
    }
    
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2]) != numVert*2 ) {
        mexErrMsgIdAndTxt("ARP:mRigidBodyProperties:size","p must be size numvertx2.");
    }
    if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetNumberOfElements(prhs[3]) != numVert*2 ) {
        mexErrMsgIdAndTxt("ARP:mRigidBodyProperties:size","pdot must be size numvertx2.");
    }
    if( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mxGetNumberOfElements(prhs[4]) != numVert*2 ) {
        mexErrMsgIdAndTxt("ARP:mRigidBodyProperties:size","lumped mass diag must be size numvertx2.");
    }

    int numRigid = (int) mxGetScalar( prhs[0] );
    // note that if these are made in the other mex file, an int32 type seems likely here
    int* rigidIDbyVert = (int*)mxGetData( prhs[1] );
    // these come directly from the mesh and will be double
    double *particlePos = mxGetDoubles( prhs[2] );
    double *particleVel = mxGetDoubles( prhs[3] );
    double *diagMass = mxGetDoubles( prhs[4] );
    
    // Allocate space for the rigid bodies out
    plhs[0] = mxCreateDoubleMatrix( 2, numRigid, mxREAL );
    double *com = mxGetDoubles(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix( 2, numRigid, mxREAL );
    double *comdot = mxGetDoubles(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix( 1, numRigid, mxREAL );
    double *mass = mxGetDoubles(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix( 1, numRigid, mxREAL );
    double *rotMass = mxGetDoubles(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix( 1, numRigid, mxREAL );
    double *omega = mxGetDoubles(plhs[4]);
    plhs[5] = mxCreateDoubleMatrix( 2*numVert, 1, mxREAL );
    double *VertexDisp = mxGetDoubles(plhs[5]);

    //initial values to prevent NaNs
    for (int i = 0; i < numRigid; ++i) {
        com[i*2] = 0;
        com[i*2+1] = 0;
        comdot[i*2] = 0;
        comdot[i * 2+1] = 0;
        omega[i] = 0;
        rotMass[i] = 0;
        mass[i] = 0;
    }
    
    for ( int i = 0; i < numVert; ++i) {
        int rid = rigidIDbyVert[i] - 1; // convert 1 index to zero index
        if ( rid < 0 ) continue;
        double m = diagMass[i*2];
        mass[rid] += m;
        com[rid*2] += particlePos[i*2] * m;
        com[rid*2+1] += particlePos[i*2+1] * m;
        comdot[rid*2] += particleVel[i*2] * m;
        comdot[rid*2+1] += particleVel[i*2+1] * m;
    }
    // did mass weighted sum, divide by total to get body pos and vel
    for ( int rid = 0; rid < numRigid; ++rid) {
        com[rid*2] /= mass[rid];
        com[rid*2+1] /= mass[rid];
        comdot[rid*2] /= mass[rid];
        comdot[rid*2+1] /= mass[rid];
    }
    // easier to pick up display positions, J and omega in a second pass
    // compute omega so as to conserve angular momentum
    for ( int i = 0; i < numVert; ++i) {
        int rid = rigidIDbyVert[i] - 1; // convert 1 index to zero index
        if ( rid < 0 ) continue;
        double m = diagMass[i*2];
        double x = particlePos[i*2] - com[rid*2];
        double y = particlePos[i*2+1] - com[rid*2+1];
        VertexDisp[i*2] = x;
        VertexDisp[i*2+1] = y;
        rotMass[rid] += m*(x*x+y*y);
        double dxvy = x * (particleVel[i*2+1] - comdot[rid*2+1]);
        double dyvx = y * (particleVel[i*2] - comdot[rid*2]);        
        omega[rid] += m * (dxvy - dyvx); // cute 2x2 cross product trick        
    }
    // a final pass to fix up our mass weighted computation of omega (i.e., preserving angular momentum)
    for ( int rid = 0; rid < numRigid; ++rid) {
        omega[rid] /= rotMass[rid];
    }
}
