
%[Vgm,Tgm,Fgm] = readMSH('icosasphere.msh'); % gmsh default save?
%[Vgm,Tgm,Fgm] = readMESH('icosasphere-m.mesh'); % gmsh export inria medit format?
% icosaspherematlab % export for matlab
% Vgm = msh.POS;
% Tgm = msh.TETS(:,1:4);


% [Vtw,Ttw,Ftw] = readMESH('icosasphere.mesh'); % tetwild 
% 
% V = readNODE( 'data3d/icosasphere.1.node' );
% T = readELE( 'data3d/icosasphere.1.ele' );
% 
% T = Ttg;
% V = Vtg;


% [V,T] = readMESH('capsule.mesh'); % tetwild 


V = readNODE( 'cubeicosadiff.1.node' );
T = readELE( 'cubeicosadiff.1.ele' );


MakeVideo = 1;
   
if ( MakeVideo==1 )
    v = VideoWriter("out/" + "cutaway-cid1",'MPEG-4');
    open(v);
end

c = repmat([0 1 0], size(T,1),1 );

camlight('right', 'infinite');
%tm = tetramesh( T,V,c, 'FaceAlpha', 0.6 );
tm = tetramesh( T,V, 'FaceAlpha', 0.9 );
camlight('right', 'infinite');
ylabel('y')
xlabel('x')
zlabel('z')

xmin = min(V(:,1));
xmax = max(V(:,1));
ymin = min(V(:,2));
ymax = max(V(:,2));
zmin = min(V(:,3));
zmax = max(V(:,3));

axis( [xmin xmax ymin ymax zmin zmax] );

xcut = linspace(xmin,xmax,60);

for x = xcut

    for i = 1:numel(tm)
        if ( any( tm(i).Vertices(:,2) < x ) )
            tm(i).Visible = 'off';
        end
    end
    if MakeVideo == 1
        F = getframe(gcf);
        writeVideo(v,F);
    else
        drawnow;
    end
end
for x = fliplr(xcut)

    for i = 1:numel(tm)
        if ( any( tm(i).Vertices(:,2) > x ) )
            tm(i).Visible = 'on';
        end
    end
    if MakeVideo == 1
        F = getframe(gcf);
        writeVideo(v,F);
    else
        drawnow;
    end
end

if MakeVideo == 1
    close(v);
end


