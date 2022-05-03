function [W, V2Htet] = meshSkinning(Vsurface,TetMesh3D)
    % Map each V to a tet(H_tet) and corresponding barycentric coordinates
    V = reshape(TetMesh3D.p',3,[])';
    V2Htet = zeros(1, size(Vsurface,1));
    V2BaryCentric = zeros(size(Vsurface,1), 4);
    
%     debugMinDist = zeros(size(Vsurface,1),1);
    for ii=1:size(Vsurface,1)
        minDistance = Inf;
        localSurfaceV = Vsurface(ii,:);
        for jj=1:size(TetMesh3D.t,1)
            % candidate tets vertices
            tetIds = TetMesh3D.t(jj,:);
            cand_TV = V(tetIds',:);     % 4x3
            % try to compute the barycentric coordinate of V(ii,:) in cand_TV
            centroid = mean(cand_TV,1);
            distanceVector = localSurfaceV - centroid;
            dist = norm(distanceVector);
            if dist < minDistance %assigns the vertex to the closest Tet
                minDistance = dist;
%                 debugMinDis(ii) = dist;
                V2Htet(ii) = jj;
                bcCoeff = barycentric_coordinates(localSurfaceV,cand_TV(1,:),cand_TV(2,:),cand_TV(3,:),cand_TV(4,:));       
                V2BaryCentric(ii,:) = bcCoeff;
            end
        end

    end

%     scatter3(V(:,1),V(:,2),V(:,3));
%     hold on;
%     scatter3(Vsurface(:,1),Vsurface(:,2),Vsurface(:,3),'MarkerFaceColor',[0 .75 .75]);
    
    assert(all(V2Htet~=0));
    
    % build barycentric J matrix such that BJ * SV(:,k) = V(:,k)
    row = size(Vsurface,1);
    col = size(V,1);
    W = sparse(row, col);
    for ii=1:size(Vsurface,1)
        W(ii,TetMesh3D.t(V2Htet(ii),:)) = V2BaryCentric(ii,:);
    end
end

