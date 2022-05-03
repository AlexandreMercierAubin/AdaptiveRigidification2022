function saveSVG(filename,V,E,F,S,W)
  %
  % saveSVG(filename,V,E,F,S,W)
  %
  % Inputs:
  %   V  #V by 2 list of vertex positions
  %   E  #E by 2 list of edge indices into V
  %   F  #3 list of decimal fill color values
  %   S  #3 list of decimal stroke color values
  %   W  stroke width values
  
  P = [];
  for e = E'
    P = [P ; V(e(2),:) - V(e(1),:)];
  end
  
  fh = fopen(filename,'w');
  fprintf(fh,['<?xml version="1.0" encoding="UTF-8" standalone="no"?> <svg version="1.1" id="Layer_1" ', ...
    ' xmlns="http://www.w3.org/2000/svg" ', ...
    ' xmlns:xlink="http://www.w3.org/1999/xlink" ', ...
    ' width="%d" height="%d" ', ...
    ' xml:space="preserve">\n'], max(V)-min(V));
  fprintf(fh, ...
    ['<path d="m']);
  fprintf(fh, ...
    [' %f,%f '], ...
    P');
  fprintf(fh, ...
    ['z" style="fill:rgb(%d,%d,%d);stroke:rgb(%d,%d,%d);stroke-width:%d" />'], ...
    [round(255*[F S]) W]');
  fprintf(fh,'</svg>\n');
  fclose(fh);
end