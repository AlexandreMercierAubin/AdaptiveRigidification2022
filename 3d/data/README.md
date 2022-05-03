# Using and Making Tetrahedral Meshes

This directory contains a nubmer of 3D triangular and tetrahedral meshes.
The puprpose of this document is to give a general guide for successful 
mesh creation.

## Starting with good meshes

While tetwild can accomodate bad meshes, and gmsh can be set up (perhaps 
with some difficulty) to reparameterize the sruface before meshing), it 
is often good to start with a high quality surface mesh.

- well shaped triangles (all close to equilateral)
- triangles of equal size 
- single closed manifold

When building examples from scratch, there is the opportunity to take these
considerations into account.  When working with a bad mesh, it might be 
good advice to give up, look elsewhere, or start over, rather than trying
to repair.  Or otherwise, consider tetwild with embedding weights as a 
default solution?

- Blender
- openscad
- http://graphite.gforge.inria.fr/  
- meshlab
- openflipper

## Compiling on windows

The bare minimum is visual studio community and cmake.

- https://visualstudio.microsoft.com/vs/community/
- https://docs.microsoft.com/en-us/cpp/build/cmake-projects-in-visual-studio?view=msvc-160
- Alternatively https://cmake.org/

## Compiling Tetgen on windows

After downloading and unzipping code, the standard process on windows
is almost identical to the steps on linux or OSX.  The only trick is
to remember to use msbuild instead of make, and to request a Release
build rather than Debug.

- Start "Developer Command Prompt for VS 2019" or equivalent
- mkdir build
- cd build
- cmake ..
- msbuild Project.sln /property:Configuration=Release

## Setting path for convenience in creating examples

Use addMeshersToWindowsPath.bat to run commands or scripts on the command
line to regenerate different tet meshes (tetgen, 
tetwild, gmsh) from surface meshes (obj, stl).

## Examples


### TetWild

TetWild --input icosasphere.stl --output icosasphere.mesh

### Tetgen (1.6.0)

tetgen -p icosasphere.stl

Or preserve quality, and impose a max volume...

tetgen -pq1.414a.001 icosasphere.stl

Or with a PLC with different materials and holes, consider an stl created
in blender.  For instance...

[V,F,N] = readSTL('sandwich.stl');
H = []; % no holes
writePOLY_tetgen('sandwich.poly',V,F,H);
 
Then, edit the file to add regions into part 4:

# Part 4 - region list
3
1  0 -0.5  0  1 
2  0  0.0  0  2
3  0  0.5  0  3

Then create the tet mesh, and note the "A" parameter is necessary to make 
sure that attributes are assigned to the tets:

tetgen -pq1.414a.01A sandwich.poly

Then check it out:

[E,A] = readELE('sandwich.1.ele');
V = readNODE('sandwich.1.node');
tetramesh( E, V, A );



### GMSH

How to set up a .geo file for the basic case?  Is there a reader for the 
default output file?  Only thing that seems easy right now is matlab export.

##