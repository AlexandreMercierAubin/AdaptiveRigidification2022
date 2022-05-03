/* WIP test code for point triangle continuous collision detection
 *========================================================*/

#include "mex.h"
#include "blas.h"
// #include <ccd.hpp>
#include <minimum_separation_root_finder.hpp>

/* 
 * The gateway function
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    
    double *ipt = mxGetDoubles(prhs[0]);
    const Eigen::Vector3d pt0(ipt[0],ipt[1],ipt[2]);
    
    double *ipt2 = mxGetDoubles(prhs[1]);
    const Eigen::Vector3d t0t0(ipt2[0],ipt2[1],ipt2[2]);
    
    double *ipt3 = mxGetDoubles(prhs[2]);
    const Eigen::Vector3d t1t0(ipt3[0],ipt3[1],ipt3[2]);
    
    double *ipt4 = mxGetDoubles(prhs[3]);
    const Eigen::Vector3d t2t0(ipt4[0],ipt4[1],ipt4[2]);
    
    double *ipt5 = mxGetDoubles(prhs[4]);
    const Eigen::Vector3d pt1(ipt5[0],ipt5[1],ipt5[2]);
    
    double *ipt6 = mxGetDoubles(prhs[5]);
    const Eigen::Vector3d t0t1(ipt6[0],ipt6[1],ipt6[2]);
    
    double *ipt7 = mxGetDoubles(prhs[6]);
    const Eigen::Vector3d t1t1(ipt7[0],ipt7[1],ipt7[2]);
    
    double *ipt8 = mxGetDoubles(prhs[7]);
    const Eigen::Vector3d t2t1(ipt8[0],ipt8[1],ipt8[2]);
    
    double dist = mxGetScalar(prhs[8]);
    
    double toi = 0.000001;
//     bool isHit = ipc::point_triangle_ccd_broadphase(pt0, t0t0, t1t0, t2t0, pt1, t0t1, t1t1, t2t1, dist);
//     bool isHit = ipc::point_triangle_ccd(pt0, t0t0, t1t0, t2t0, pt1, t0t1, t1t1, t2t1, dist);
//     bool isHit = eccd::vertexFaceCCD(pt0, t0t0, t1t0, t2t0, pt1, t0t1, t1t1, t2t1);
//     bool isHit = ccd::vertexFaceCCD(pt0, t0t0, t1t0, t2t0, pt1, t0t1, t1t1, t2t1, ccd::BSC);
    bool isHit = msccd::root_finder::vertexFaceMSCCD(pt0, t0t0, t1t0, t2t0, pt1, t0t1, t1t1, t2t1, dist,toi);
    
    // create the output vectors 
    plhs[0] = mxCreateLogicalScalar(isHit);
}
