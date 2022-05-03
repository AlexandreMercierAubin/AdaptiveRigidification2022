% Compile the cpp mex code needed for 2D and 3D adaptive rigidificaiton
rootPwd = pwd;
cd 2d/mex
mex -R2018a mexComputeSTVKGradHess2D.cpp
mex -R2018a mexEdotNorm2D.cpp
mex -R2018a mexEdiffNorm2D.cpp
mex -R2018a mexFindMatches2D.cpp
mex -R2018a mexPGS2D.cpp
mex -R2018a mexPGS2DwithJAinvJT.cpp
mex -R2018a mexRigidBodyProperties2D.cpp
mex -R2018a mexRigidConnectedComponents2D.cpp
cd ../..
cd 3d/mex
mex -R2018a mexSTVK3D.cpp
mex -R2018a mexNeoHookean3D.cpp
mex -R2018a mexCorotational3D.cpp
mex -R2018a mexEdotNorm3D.cpp
mex -R2018a mexFindMatches3D.cpp
mex -R2018a mexPGS3D.cpp
mex -R2018a mexPGS3DwithJAinvJT.cpp
mex -R2018a mexRigidBodyProperties3D.cpp
mex -R2018a mexRigidBodyConnectedComponents3D.cpp
cd ../..
cd util/mex
mex -R2018a mMD5.c
cd ../..
%-------you can temporarily comment this if you haven't built ipc-toolkit
%ipc or ECCD won't work
% ipctoolkitIncludePath = ['-I' fullfile(rootPwd, 'lib','ipc-toolkit','build','include')];
% ipctoolkitLibPath = ['-L' fullfile(rootPwd, 'lib','ipc-toolkit','build','Release')];
% ipctoolkitlibPath = ['-l','IPCToolkit.lib'];
% eigenpath = ['-I' fullfile(rootPwd, 'lib','eigen')];
% mex('-v', '-R2018a', eigenpath, ipctoolkitIncludePath, 'pointTriangleCCD.cpp', ipctoolkitLibPath, ipctoolkitlibPath);
% mex('-v', '-R2018a', eigenpath, ipctoolkitIncludePath, 'mexIPC.cpp', ipctoolkitLibPath, ipctoolkitlibPath);

%-------

%cd ../..
