# Adaptive Rigidification project code for Siggraph 2022
Matlab code for the adaptive rigidification project. This code base provides an implementation of the adaptive rigidification technique as seen at Siggraph 2022. Matlab is intuitive which makes the project easy to setup and results quickly reproducable. To test the various scenes featured in the paper, simply follow the instructions below then open an example file and click play in the editor. Alternatively the name of the example can be called in the command prompt which will launch the simulation. The scenes feature heterogenous materials, multiple objects, scripted animations, surface mapping and more!

## How to use

### Dependencies
required: opengl, visual studio (or any C++ compilers), cmake, matlab, gptoolbox(for 3D)
GPtoolbox might require: 
	embree https://github.com/embree/embree
	libigl https://github.com/libigl/libigl
	boost 1.48 https://www.boost.org/users/history/version_1_48_0.html
	to be detected automatically on windows use the path C:\local\boost_1_48_0\boost\*
Tip: install visual studio before cmake so its compiler is automatically added to cmake

GpToolBox is already included in the repo as a subrepo. Simply call `git submodule update --init --recursive`.
cmake is needed to build gptoolbox's mex. You will need the boost library. 
I recommend not to build gptoolbox with eltopo if you are on windows.

We depend on GPToolbox for contacts(signed distances) as well as loading some mesh files.

Some matlab addons are needed or recommended such as Deep Learning Toolbox, Computer Vision Toolbox, Matlab Compiler, Signal Processing Toolbox, and Image processing Toolbox. Those can be found and downloaded from matlab's add-on menu.

### Installation
#### run the following commands in a command prompt
```
git submodule update --init --recursive
#if the recursive update fails then manually go to lib/gptoolbox/mex and build
cd lib/gptoolbox/mex
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target install --config Release
```

#### run the following commands in the matlab console
Once the project is set up, add the files to path by calling `addAdaptiveRigidificationToPath();` in the command prompt.
Followed by `mexAdaptiveRigidification();` to build the project's mex code.
```
addAdaptiveRigidificationToPath();
savepath();
mexAdaptiveRigidification();
```
By now you should be all set.

## Running scenes
Try running examples in 2d/examples or 3d/examples.
For instance while being in the root folder, you could open 2d/examples/cantilever_pinned.m and click play or simply call cantilever_pinned in the command prompt.
Why don't you try creating your own scene and play with it?

### FAQ
If you get a signed-distance error when running the 3D code, then you probably did not set up the compiled mex from gptoolbox properly.

Feel free to contact me if you find bugs and to send pull requests with fixes. 
We likewise have a private repository that we maintain with new research. 
Contact us to gain access to the private repository if you wish to collaborate with us on our newest/most exciting stuff.
Please cite us when using this simulator in your own research.

Special thanks to the collaborators Alexander Winter, David IW Levin and Paul Kry who participated in creating this simulator.
This simulator is now used by many students of the McGill physics-based animation lab .
We wish to keep expanding it over time.

