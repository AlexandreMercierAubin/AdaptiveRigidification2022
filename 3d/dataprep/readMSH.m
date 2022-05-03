function [V,T,F] = readMSH(filename)
  % Input:
  %   filename path to gmsh .msh file
  % Outputs:
  %   V  #V by 3 list of vertex positions
  %   T  #T by 4 list of tet indices
  %   F  #F by 3 list of triangle indices
   
  function expect_line(f,str)
   line = fscanf(f,'%s\n',1);
   if strcmp(line,str) ~= 1
     assert(false,'"%s" ~= "%s"',line,str);
   end
  end

  f = fopen(filename,'r');
  expect_line(f,'$MeshFormat');
  vfd = fscanf(f,'%g %g %g\n',3);
  version = vfd(1);
  filetype = vfd(2);
  datasize = vfd(3);
  assert(datasize == 8, 'only doubles are supported');
  ascii = filetype == 0;
  if ascii
   expect_line(f,'$EndMeshFormat');
   expect_line(f,'$Nodes');
   
   n = fscanf(f,'%d\n',1);
   V = zeros(1,3);
   for i = [1:n]
    id = fscanf(f,'%d\n',1);
    V(id,1) = fscanf(f,'%f\n',1);
    V(id,2) = fscanf(f,'%f\n',1);
    V(id,3) = fscanf(f,'%f\n',1);
   end
   expect_line(f,'$EndNodes');
   expect_line(f,'$Elements');
   n = fscanf(f,'%d\n',1);
   T = zeros(1,4);
   F = zeros(1,3);
   nbTet=0;
   for i = [1:n]
    id = fscanf(f,'%d\n',1);
    size = fscanf(f,'%d\n',1);
    tmp = fscanf(f,'%d %d %d\n',3);
    if (size==4)
        T(id,1) = fscanf(f,'%f\n',1);
        T(id,2) = fscanf(f,'%f\n',1);
        T(id,3) = fscanf(f,'%f\n',1);
        T(id,4) = fscanf(f,'%f\n',1);
        nbTet=i;
    else 
        F(id-nbTet,1) = fscanf(f,'%f\n',1);
        F(id-nbTet,2) = fscanf(f,'%f\n',1);
        F(id-nbTet,3) = fscanf(f,'%f\n',1);
    end
   end
   expect_line(f,'$EndElements');
  else
   assert(false,"Not implemented.")
  end
  fclose(f);
end
