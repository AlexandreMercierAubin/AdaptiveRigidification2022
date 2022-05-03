To render the animation, we first visualize them in the matlab editor. 
Once we achieve visually interesting motions, we export per frame meshes in either obj or ply format.
We then use stop-motion-obj to import the frames into blender and use it to give that extra visual quality to the animations. 
For larger scenes, I recommend using the streaming import mode and then render separate png files.
You should execute the render from commandline as blender tends to crash when doing it with stop-motion-obj in the application.
Once you have a bunch of images, simply send them to ffmpeg and there you go, you have a mp4 ready to be showed.

In 2d meshes are generated using triangle and 3d mostly use tetwild.

## Using FFMPEG reminder

If rendering to png, this can convert to avi
ffmpeg -i %4d.png -preset veryslow -qscale 0 -vcodec libx264 -pix_fmt yuv420p -f mp4 out.mp4

-i is the filename %nd means that we have trailing numbers of length n at the location specified ex: eman%3dname.png will start with eman001name.png
-vcodec should always be define otherwise it will default to a poor quality mpeg-1. You probably want either mpeg4 or libx264
-f file extension

## Tetwild command reminder

./TetWild --input [input] --output [output] -l [ideal edge-length ex:0.1]