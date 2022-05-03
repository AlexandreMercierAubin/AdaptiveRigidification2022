function exportobj(p, t, fname)

fileID = fopen( fname, 'w' );
for i = 1:size(p,1),
    fprintf( fileID, 'v %f %f 0\n', p(i,1), p(i,2) );
end;
for i = 1:size(t,1),
    fprintf( fileID, 'f %d %d %d\n', t(i,1), t(i,2), t(i,3) );
end;
fclose(fileID);