/*==========================================================
 * mRigidConnectedComponents.c
 *
 * Computes connected rigid components based on Edot history, but likewise 
 * prevents rigidification in any elements that are currently showing a 
 * large approxEdot from the quick solve.  This method uses the triangle
 * connectivity graph (possibly just edge adjacency, but could also be one
 * of two other alternatives to avoid hinges).
 *
 * The general strategy does a BFS starting form the towPinTris elements, 
 * but then simply walking over all tris afterward and restarting the BFS
 * for a rigid component at any elements that have not yet been visited.
 * Temp memory is needed to keep track of vertices and triangles visited.
 * Rigidificaiton rules will combine the input EdotHistory, EdotApprox, 
 * along with one of the following:
 * 1) Tri-adj across edges, and a first come first serve greedy approach to
 * giving a vertex to a rigid component. 
 * 2) Tri-adj across edges and vertices, recognizing that hinges will be 
 * welded
 * 3) a more complex multi pass approach of identiyfing vertices that have 
 * all triangles ready to rigidfy, then allowing triangles that have all 
 * verts ready to go become rigid (this way tri-adj across edges and 
 * vertices will weld, but all surrounding triangles will be in a low Edot
 * state so it will be 100% Kosher!
 * 
 * [ numRigid, rigidIDbyTri rigidIDbyVert, isElasticVert, isBoundaryVert, isElasticTri ] = mRigidConnectedComponents( ... 
 *     EdotHistory, EdotApprox, adjMatrix, triVertIDint32, numVerts, pinnedVerts ) 
 *
 *  inputs:
 *
 *      EdotHistory     #Tri by 3? Edot history of elements
 *      EdotApprox      #tri by 1 result from quick solve no constraints
 *      adjMatrix       #tri by #tri SPARSE symmetric matrix of triangle adj
 *      triVertIDint32  #tri by 3 (4 for tets) indices (starting from 1) of verts of each element
 *      numVerts        scalar number of vertices
 *      ispinnedVerts     #vert by 1 logical identifying pinned vertices
 *      verticeValance  #vert by 1 array of the valances
 *      isStableVerts   #vert by 1 logical identifying pinned vertices
 *
 *  outputs:
 *
 *      numRigid        scalar, number of rigids (possibly zero)
 *      rigidIDbyTri	#Tri by 1 vector of rigid IDs 
 *      rigidIDbyVert   #Vert by 1 vector of rigid IDs
 *      isVertElastic   #Vert by 1 logical identifying rigid verts (i.e., rigidIDbyVert == -1, and avoid a find on the matlab side)
 *      isVertBoundary  #Vert by 1 logical identifying rigid verts on the boundary of a rigid
 *      isTriElastic    #Tri by 1 logical identifying rigid tris (i.e., rigidIDbyTri == -1, and avoid a find on the matlab side) 
 *
 *
 * To compile type: mex -R2018a mexRigidConnectedComponents2D.cpp
 *========================================================*/

#include "mex.h"
#include "blas.h"
#include <map>
#include <vector>
#include <queue>

void visitComponent(const size_t rootIndex, mxLogical* isVertBoundary, mxLogical* isTriElastic,
                    mxLogical* isVertElastic,int* verticesRigidCount, int* rigidIDbyVert,
                    int* rigidIDbyTri, mwIndex* jc, mwIndex* ir, double* adjMatrix, int* triVertIDint32,
                    const size_t numTriangles, bool* pinnedVerts, 
                    int &numRigid, const size_t vSize, std::queue<size_t> &toVisit, std::vector<bool> &isComponentPinned, bool* isStableTri){


       if(rigidIDbyTri[rootIndex] == 0){
        numRigid++;
        isComponentPinned.push_back(false);
        toVisit.push(rootIndex);
        while(!toVisit.empty()){ 
            size_t tId = toVisit.front();
            toVisit.pop();

            //check for hinges
            bool hasHinge = false;
            int pinnedCount = 0;
            for (size_t v = 0; v < vSize; ++v) {
                int vertex = triVertIDint32[tId + v * numTriangles] - 1;

                if (pinnedVerts[vertex]) {
                    pinnedCount++;
                }

                //The vertex was already claimed by a different rigid body or must remain elastic, exit. Same for a single pinned dof not tied to a stable tri
                bool isClaimed = (rigidIDbyVert[vertex] > 0 && rigidIDbyVert[vertex] != numRigid);
                bool isNotStable = pinnedCount > 0 && !isStableTri[tId] && !isComponentPinned[numRigid-1];
                if (isClaimed || isVertElastic[vertex] || isNotStable) {
                    hasHinge = true;
                    break;
                }
            }
            
            if (hasHinge) {
                rigidIDbyTri[tId] = -1;
                isTriElastic[tId] = true;
                
                // if this is the first triangle explored in a new rigid, then delete the rigid body
                if (tId == rootIndex) {
                    numRigid--;
                    isComponentPinned.pop_back();
                }
                continue;
            }

            rigidIDbyTri[tId] = numRigid;
            isTriElastic[tId] = false;
            size_t starting_row_index = jc[tId];
            size_t stopping_row_index = jc[tId + 1];

            //No hinge, nice! claim the vertices!
            
            for(size_t v = 0; v < vSize; ++v){
                int vertex = triVertIDint32[tId + v * numTriangles]-1; // index starts at 0 instead of 1
                rigidIDbyVert[vertex] = numRigid;
                verticesRigidCount[vertex]++;
            }
            if(pinnedCount>0){
                isComponentPinned[numRigid - 1] = true;
            }

            for ( size_t current_row_index = starting_row_index; current_row_index < stopping_row_index; current_row_index++ )  {
                size_t row = ir[current_row_index];
                
                if(rigidIDbyTri[row] == 0){
                    rigidIDbyTri[row] = numRigid;
                    toVisit.push(row);
                }
            }
        }
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    
        /* check for proper number of arguments */
    if ( nrhs != 9 ) {
        mexErrMsgIdAndTxt("ARP:rigidConnectedComponents:nrhs","9 inputs required.");
    }
    if ( nlhs != 7 ) {
        mexErrMsgIdAndTxt("ARP:rigidConnectedComponents:nlhs","7 outputs required.");
    }
    /* TODO: more input checking!!! or play dangerous and assume perfect input? or is this stuff slow?*/
    bool* EdotHistory   = mxGetLogicals(prhs[0]);     //#Tri by 3? Edot history of elements, logical array
    const size_t framesObserved = mxGetN(prhs[0]);
    const size_t numTriangles = mxGetM(prhs[0]);
    bool isLogicalEdotHistory = mxIsLogical(prhs[0]);
    if (!isLogicalEdotHistory) {
        mexErrMsgIdAndTxt("ARP:rigidConnectedComponents:type", "EdotHistory must be logical");
    }

    bool* EdotApprox    = mxGetLogicals(prhs[1]);    //#tri by 1 result from quick solve no constraints , logical array
    bool isLogicalEdotApprox = mxIsLogical(prhs[1]);
    const size_t rowsEdotApprox = mxGetM(prhs[1]);
    const size_t colsEdotApprox = mxGetN(prhs[1]);
    if (!isLogicalEdotApprox) {
        mexErrMsgIdAndTxt("ARP:rigidConnectedComponents:type", "EdotApprox must be logical");
    }
    if (numTriangles != rowsEdotApprox) {
        mexErrMsgIdAndTxt("ARP:rigidConnectedComponents:size", "EdotApprox: Expecting a column vector of size Tri");
    }

    double* adjMatrix   = mxGetDoubles(prhs[2]);    //#tri by #tri SPARSE symmetric matrix of triangle adj (3xentries)
    mwIndex* ir = mxGetIr( prhs[2] );
    mwIndex* jc = mxGetJc( prhs[2] );
    if (numTriangles != mxGetM(prhs[2]) || numTriangles != mxGetN(prhs[2])) {
        mexErrMsgIdAndTxt("ARP:rigidConnectedComponents:size", "adjMatrix: Expecting a matrix Tri by Tri");
    }
    if (!mxIsDouble(prhs[2])) {
        mexErrMsgIdAndTxt("ARP:rigidConnectedComponents:type", "adjMatrix: Expecting matrix to be doubles");
    }
    
 	int* triVertIDint32 = (int*)mxGetData(prhs[3]);       //#tri by 3 indices (starting from 1) of verts of each element
    const size_t vSize = mxGetN(prhs[3]);
    if (numTriangles != mxGetM(prhs[3]) || 3 != vSize) {
        mexErrMsgIdAndTxt("ARP:rigidConnectedComponents:size", "triVertIDint32: Expecting a matrix Tri by 3");
    }
    
 	double numVerts = mxGetScalar(prhs[4]); //scalar number of vertices
    if (!mxIsScalar(prhs[4])) {
        mexErrMsgIdAndTxt("ARP:rigidConnectedComponents:type", "numVerts must be scalar");
    }

    bool* pinnedVerts = mxGetLogicals(prhs[5]);    //#vert by 1 logical identifying pinned vertices
    const size_t rowsPinnedVerts = mxGetM(prhs[5]);
    bool isLogicalPinned = mxIsLogical(prhs[5]);
    if (!isLogicalPinned) {
        mexErrMsgIdAndTxt("ARP:rigidConnectedComponents:type", "PinnedVerts must be logical");
    }
    if (rowsPinnedVerts != numVerts) {
        mexErrMsgIdAndTxt("ARP:rigidConnectedComponents:size", "PinnedVerts: Expecting a column vector of size numVerts");
    }

    double* verticeValance = mxGetDoubles(prhs[6]);    //#vert by 1 valence of vertices
    if (!mxIsDouble(prhs[6]) || mxGetM(prhs[6]) != numVerts){
        mexErrMsgIdAndTxt("ARP:rigidConnectedComponents:size", "verticeValance must be size numvert.");
    }
    
    double* stablePinnedTris = mxGetDoubles(prhs[7]);
    const size_t numStableTris = mxGetM(prhs[7]);
    const size_t colsStableTris = mxGetN(prhs[7]);
    if (numStableTris < colsStableTris && numStableTris != 0) {
        mexErrMsgIdAndTxt("ARP:rigidConnectedComponents:size", "stablePinnedTris: Expecting a column vector");
    }
    if (!mxIsDouble(prhs[7])) {
        mexErrMsgIdAndTxt("ARP:rigidConnectedComponents:type", "stablePinnedTris: Expecting array to be doubles");
    }

    bool* isStableTri = mxGetLogicals(prhs[8]);
    bool isLogicalStableTri = mxIsLogical(prhs[8]);
    const size_t rowsisStableTri = mxGetM(prhs[8]);
    if (!isLogicalStableTri) {
        mexErrMsgIdAndTxt("ARP:rigidConnectedComponents:type", "isStableTri must be logical");
    }
    if (rowsisStableTri != numTriangles) {
        mexErrMsgIdAndTxt("ARP:rigidConnectedComponents:size", "isStableTri: Expecting a column vector of size numTriangles");
    }

    //new vars for output
    plhs[1] = mxCreateNumericMatrix((mwSize)numTriangles, 1, mxINT32_CLASS, mxREAL);//rigidIDbyTri: #Tri by 1 vector of rigid IDs
    int* rigidIDbyTri = (int*)mxGetData(plhs[1]);
    plhs[2] = mxCreateNumericMatrix((mwSize)numVerts, 1, mxINT32_CLASS, mxREAL);    //rigidIDbyVert: #Vert by 1 vector of rigid IDs
    int* rigidIDbyVert = (int*)mxGetData(plhs[2]);
    plhs[3] = mxCreateLogicalMatrix((mwSize)numVerts,1);    //isVertElastic: #Vert by 1 logical identifying rigid verts (i.e., rigidIDbyVert == -1, and avoid a find on the matlab side)
    mxLogical* isVertElastic = mxGetLogicals(plhs[3]);
    plhs[4] = mxCreateLogicalMatrix((mwSize)numVerts,1);    //isVertBoundary: #Vert by 1 logical identifying rigid verts on the boundary of a rigid
    mxLogical* isVertBoundary = mxGetLogicals(plhs[4]);
    plhs[5] = mxCreateLogicalMatrix((mwSize)numTriangles,1);//isTriElastic: #Tri by 1 logical identifying rigid tris (i.e., rigidIDbyTri == -1, and avoid a find on the matlab side) 
    mxLogical* isTriElastic = mxGetLogicals(plhs[5]);  
    
    int* verticesRigidCount = new int[static_cast<int>(numVerts)]{0};
    //init elements
    for(size_t i = 0; i < numVerts; ++i){
        isVertElastic[i] = true;
    }
    
    for ( size_t i = 0; i < numTriangles; ++i ) {
        isTriElastic[i] = true;
        rigidIDbyTri[i] = -1;
        bool eligible = !EdotApprox[i];
        for(size_t frame = 0; frame < framesObserved; ++frame){
            eligible = eligible && EdotHistory[i + frame*numTriangles];
        }
        for (size_t v = 0; v < vSize; ++v) { //init vert as elastic 
            int vertex = triVertIDint32[i + v * numTriangles] - 1; // index starts at 0 instead of 1
            isVertBoundary[vertex] = eligible;
        }
        if (eligible) { 
            rigidIDbyTri[i] = 0; 

            for (size_t v = 0; v < vSize; ++v) { //set eligible vertices as rigid
                int vertex = triVertIDint32[i + v * numTriangles] - 1; // index starts at 0 instead of 1
                isVertElastic[vertex] = false;
            }
        }
    }
    
    
    for ( size_t i = 0; i < numTriangles; ++i ) {
        //Forces the elastification of all dofs of an element tagged for elastification.
        if (EdotApprox[i]){ //force elastification of all of the vertices of a triangle to elastify
            for (size_t v = 0; v < vSize; ++v) { //init vert as elastic 
                int vertex = triVertIDint32[i + v * numTriangles] - 1; // index starts at 0 instead of 1
                isVertElastic[vertex] = true;
            }
        }
    }
    
    //rigidIDbyTri acts as isRigid and isTriElastic acts as visited
    // it saves memory
    int numRigid = 0;
    std::queue<size_t> toVisit;
    std::vector<bool> isComponentPinned;
    
    //start the search with stable pinned tris
    for ( size_t i = 0; i < numStableTris; ++i ) {
        size_t t = (size_t)stablePinnedTris[i] - 1;
        visitComponent(t, isVertBoundary, isTriElastic, isVertElastic,
                        verticesRigidCount, rigidIDbyVert, rigidIDbyTri,
                        jc, ir, adjMatrix, triVertIDint32, numTriangles, 
                        pinnedVerts, numRigid, vSize, toVisit, isComponentPinned, isStableTri);
    }
    
    
    for ( size_t i = 0; i < numTriangles; ++i ) {
        visitComponent(i, isVertBoundary, isTriElastic, isVertElastic,
                        verticesRigidCount, rigidIDbyVert, rigidIDbyTri,
                        jc, ir, adjMatrix, triVertIDint32, numTriangles, 
                        pinnedVerts, numRigid, vSize, toVisit, isComponentPinned, isStableTri);
    }

    //remove elastify unclaimed vertices (removed hinge with single triangle)
    for (size_t i = 0; i < numVerts; ++i) {
        if(rigidIDbyVert[i] == 0){
            isVertElastic[i] = true;
        }
    }

    plhs[0] = mxCreateDoubleScalar(numRigid);
    
    plhs[6] = mxCreateLogicalMatrix(isComponentPinned.size(),1);//isTriElastic: #Tri by 1 logical identifying rigid tris (i.e., rigidIDbyTri == -1, and avoid a find on the matlab side) 
    mxLogical* pinnedComponent = mxGetLogicals(plhs[6]);
    for (size_t i = 0; i<isComponentPinned.size(); ++i){
        pinnedComponent[i] = isComponentPinned[i];
    }
    
    for (size_t i = 0; i<numVerts; ++i){
        if(verticeValance[i]==verticesRigidCount[i]){
            isVertBoundary[i] = false;
        }
    }
    delete[] verticesRigidCount;
}
