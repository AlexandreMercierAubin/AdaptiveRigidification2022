function [indices] = boxROI(V,F,boxes)
    % Output indices of points in V which are in boxes
    % Input:
    %   V  list of vertices
    %   boxes  aligned boxes [[minx,miny,minz,maxx,maxy,maxz],...]
    % Output:
    %   indices list of point indices in ROI
    
    indices=[];
    figure;
    plot_mesh(V,F);
    title("Fixed points");
    hold on;
    if (size(boxes,1)>1)
        for box = boxes
            indices = [indices, getIndices(V,box)];
        end
    else
        indices = [indices, getIndices(V,boxes)];
    end
end

function [indices] = getIndices(V,box)
    indices=[];
    minx=min(box(1),box(4));miny=min(box(2),box(5));minz=min(box(3),box(6));
    maxx=max(box(1),box(4));maxy=max(box(2),box(5));maxz=max(box(3),box(6));
    i=0;
    for position = V'
        i=i+1;
        if (position(1)>=minx) && (position(1)<=maxx)
            if (position(2)>=miny) && (position(2)<=maxy)
                if (position(3)>=minz) && (position(3)<=maxz)
                    indices = [indices, i];
                    scatter3(V(i,1),V(i,2),V(i,3),'r');
                end
            end
        end
    end
end
