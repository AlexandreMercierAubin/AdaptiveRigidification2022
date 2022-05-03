function saveSVG(filename,V,E,F,S,W)
  %
  % saveSVG(filename,V,E,F,S,W)
  %
  % Inputs:
  %   V  #V by 2 list of vertex positions
  %   E  #E by 2 list of edge indices into V
  %   F  #E by 3 list of decimal fill color values
  %   S  #E by 3 list of decimal stroke color values
  %   W  #E by 1 list of stroke width values
  % Example:
  %   writeSVG('wolf-rods.svg',svg.V,svg.E,svg.F,svg.S,svg.W);
  %   !cat wolf-rods.svg | sed -e "s/fill:[^;]*;//g" | /usr/local/bin/cairosvg - -o wolf-rods.pdf
  %
  % Copyright https://github.com/alecjacobson/gptoolbox/blob/master/mesh/writeSVG.m
  
  fh = fopen(filename,'w');
  fprintf(fh,['<svg version="1.1" id="Layer_1" ' ...
    ' xmlns="http://www.w3.org/2000/svg" ' ...
    ' xmlns:xlink="http://www.w3.org/1999/xlink" ' ...
    ' width="%d" height="%d" ' ...
    ' xml:space="preserve">\n'], max(V)-min(V));
  fprintf(fh, ...
    ['<path d= "m x="%f",y="%f" z' ...
    'style="fill:rgb(%d,%d,%d);stroke:rgb(%d,%d,%d);stroke-width:%d" />'], ...
    [V(E(:,1),:) V(E(:,2),:) round(255*[F S]) W]');
  fprintf(fh,'</svg>\n');
  fclose(fh);
end