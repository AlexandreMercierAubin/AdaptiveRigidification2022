/*==========================================================
 * mexRigidBodyProperties3D.cpp
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
 *      p               #vertx 3 by 1 vector of positions
 *      pdot            #vertx 3 by 1 vector of velocities
 *      mass            #vertx 3 by 1 lumped diagonal mass matrix
 *
 *  outputs:
 *
 *      com                 3 by #rigid position of center of mass of rigid
 *      comdot              3 by #rigid velocity of center of mass of rigid
 *      mass                1 by #rigid mass of a rigid
 *      rotMass             9 by #rigid rotational inertia of rigid
 *      angularMomentum     3 by #rigid current angular momentum, don't forget to solve AngularVelocity = Inertia \ angularMomentum
 *      VertexDisp          position of vertex wrt its COM
 *      
 * To compile type: mex -R2018a mexRigidBodyProperties3D.cpp
 *========================================================*/

#include "mex.h"
#include "blas.h"
#include <array>

std::array<double, 9> computemassrhatsquared(double x,double y,double z,double m){
//     double rhat[3][3]= {{0,   -z,   y},
//                           {z,    0,  -x},
//                           {-y ,  x,   0}};
    
    double mxy = m*x*y;
    double mxz = m*x*z;
    double myz = m*y*z;

    //the output will be transposed 
    return {m*(z*-z + y*-y),   mxy,              mxz,
          mxy,                 m*(z*-z+x*-x),    myz,
          mxz,                 myz,              m*(y*-y + x*-x)};
}

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
        mexErrMsgIdAndTxt("ARP:mRigidBodyProperties:size","num rigid must be a scalar.");
    }
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2]) != numVert*3 ) {
        mexErrMsgIdAndTxt("ARP:mRigidBodyProperties:size","p must be size numvertx3.");
    }
    if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxGetNumberOfElements(prhs[3]) != numVert*3 ) {
        mexErrMsgIdAndTxt("ARP:mRigidBodyProperties:size","pdot must be size numvertx3.");
    }
    if( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mxGetNumberOfElements(prhs[4]) != numVert*3 ) {
        mexErrMsgIdAndTxt("ARP:mRigidBodyProperties:size","lumped mass diag must be size numvertx3.");
    }

    int numRigid = (int) mxGetScalar( prhs[0] );
    // note that if these are made in the other mex file, an int32 type seems likely here
    int* rigidIDbyVert = (int*)mxGetData( prhs[1] );
    // these come directly from the mesh and will be double
    double *particlePos = mxGetDoubles( prhs[2] );
    double *particleVel = mxGetDoubles( prhs[3] );
    double *diagMass = mxGetDoubles( prhs[4] );
    
    // Allocate space for the rigid bodies out
    plhs[0] = mxCreateDoubleMatrix( 3, numRigid, mxREAL );
    double *com = mxGetDoubles(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix( 3, numRigid, mxREAL );
    double *comdot = mxGetDoubles(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix( 1, numRigid, mxREAL );
    double *mass = mxGetDoubles(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix( 9, numRigid, mxREAL );
    double *rotMass = mxGetDoubles(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix( 3, numRigid, mxREAL );
    double *omega = mxGetDoubles(plhs[4]);
    plhs[5] = mxCreateDoubleMatrix( 3*numVert, 1, mxREAL );
    double *VertexDisp = mxGetDoubles(plhs[5]);

    //initial values to prevent NaNs
    for (int i = 0; i < numRigid; ++i) {
        com[i * 3] = 0;
        com[i * 3 + 1] = 0;
        com[i * 3 + 2] = 0;
        comdot[i * 3] = 0;
        comdot[i * 3 + 1] = 0;
        comdot[i * 3 + 2] = 0;
        omega[i * 3] = 0;
        omega[i * 3 + 1] = 0;
        omega[i * 3 + 2] = 0;
        for(int x = 0; x<9; ++x){
            rotMass[i * 9 + x] = 0;
        }
        mass[i] = 0;
    }
    
    for ( int i = 0; i < numVert; ++i) {
        int rid = rigidIDbyVert[i] - 1; // convert 1 index to zero index
        if ( rid < 0 ) continue;
        double m = diagMass[i*3];
        mass[rid] += m;
        com[rid*3] += particlePos[i*3] * m;
        com[rid*3+1] += particlePos[i*3+1] * m;
        com[rid*3+2] += particlePos[i*3+2] * m;
        comdot[rid*3] += particleVel[i*3] * m;
        comdot[rid*3+1] += particleVel[i*3+1] * m;
        comdot[rid*3+2] += particleVel[i*3+2] * m;
    }
    // did mass weighted sum, divide by total to get body pos and vel
    for ( int rid = 0; rid < numRigid; ++rid) {
        com[rid*3] /= mass[rid];
        com[rid*3+1] /= mass[rid];
        com[rid*3+2] /= mass[rid];
        comdot[rid*3] /= mass[rid];
        comdot[rid*3+1] /= mass[rid];
        comdot[rid*3+2] /= mass[rid];
    }
    // easier to pick up display positions, J and omega in a second pass
    // compute omega so as to conserve angular momentum
    for ( int i = 0; i < numVert; ++i) {
        int rid = rigidIDbyVert[i] - 1; // convert 1 index to zero index
        if ( rid < 0 ) continue;
        double m = diagMass[i*3];
        double x = particlePos[i*3] - com[rid*3];
        double y = particlePos[i*3+1] - com[rid*3+1];
        double z = particlePos[i*3+2] - com[rid*3+2];
        VertexDisp[i*3] = x;
        VertexDisp[i*3+1] = y;
        VertexDisp[i*3+2] = z;
        
        std::array<double, 9>  inertia = computemassrhatsquared(x,y,z,m);
        for(int x = 0; x<9; ++x){
            rotMass[rid * 9 + x] -= inertia[x];
        }
        
        double vx = (particleVel[i*3]-comdot[rid*3]);
        double vy = (particleVel[i*3+1]-comdot[rid*3+1]);
        double vz = (particleVel[i*3+2]-comdot[rid*3+2]);
        
        double crossx = m*(y*vz - z*vy);
        double crossy = m*(z*vx - x*vz);
        double crossz = m*(x*vy - y*vx);
        
        omega[rid * 3] += crossx;
        omega[rid * 3+1] += crossy;
        omega[rid * 3+2] += crossz;    
    }
}
