function [W, T, SV, V, V2Htet] = voxelizer(V,FACES,side)
    % Note: # voxels is for the x direction... can end up with many more in
    % other directions!  Some models will fail without adjusting the divisions
    % and padding by a small amount.
    % [W,T,SV,V] = voxelizer(V,FACES);
    % W: weight matrix from vertices to the surface
    % T tetrahedrons of the generated tetrahedralized voxel
    % SV: surface vertices
    % V: Voxel vertices
    
    [W,BC,DV,Q] = voxelize(V,FACES,side,'Pad',1);

    % Set up a nice axis for viewing
    ll = min(DV);
    ur = max(DV);
    center = (ll+ur)/2;
    radius = (ur-ll)/2;
    buf = max(radius)*.2;
    axis( [ ...
        center(1) - radius(1) - buf, center(1) + radius(1)*4 + buf, ...
        center(2) - radius(2) - buf, center(2) + radius(2)*1 + buf, ...
        center(3) - radius(3) - buf, center(3) + radius(3)*1 + buf ] );

    % Build elements by referencing corners 
    wix = find(W==1);
    [theta,j,k] = ind2sub(size(W),wix);
    sWp1 = size(W) + [1,1,1]; % size W plus 1       % each row of el is a tuple of size 8, which is 8 indices into DV
    el = [ sub2ind( sWp1, theta  , j  , k  ), ...
           sub2ind( sWp1, theta  , j  , k+1), ...
           sub2ind( sWp1, theta  , j+1, k  ), ...
           sub2ind( sWp1, theta  , j+1, k+1), ...
           sub2ind( sWp1, theta+1, j  , k  ), ...
           sub2ind( sWp1, theta+1, j  , k+1), ...
           sub2ind( sWp1, theta+1, j+1, k  ), ...
           sub2ind( sWp1, theta+1, j+1, k+1) ];

    % Get the hexahedral mesh nodes, and hexahedral element definitions
    [SV,IM] = remove_unreferenced(DV,el);   % IM is a mapping from old index to new index
    H_hex = IM(el);  % Hexahedral elements

    % For each voxel defined by SV, H, create a regular_tet voxel and map the
    % vertices. Each voxel will be divided into 6 tets
    H_tet = zeros(size(H_hex,1)*6, 4);
    for idx=1:size(H_hex,1)
        VV = SV(H_hex(idx,:),:);        % voxel vertices
        [V_reg,T_reg,F_reg] = regular_tetrahedral_mesh(2,2,2,'Method','reflectionally-symmetric');  % a single tet-lized voxel [0,1]
        voxel_size = max(VV) - min(VV);     % should be the same in all 3 dimensions

        V_reg = V_reg .* voxel_size;        % scale and translate V_reg to VV
        V_reg = V_reg + mean(VV) - mean(V_reg);

        epsilon = mean(voxel_size) * 1e-3;
        % now try to match V_reg to VV
        Vreg2VV = zeros(1, 8);
        for ii=1:8
            Vreg2Assign = V_reg(ii,:);
            for jj=1:8
                if(norm(Vreg2Assign - VV(jj,:)) < epsilon)
                    Vreg2VV(ii) = H_hex(idx,jj);
                    break
                end
            end
        end

        % now map T_reg to indices into H
        H_tet(idx*6-5:idx*6,:) = Vreg2VV(T_reg);

    end
    
    % Next try to rebuild V from SV
    % first compute grid spacing h
    ll = min(DV);
    ur = max(DV);
    sW = size(W);
    h = (ur-ll)./sW([2,1,3]); % the hexahedron size should be the same in all dimensions!

    % Location of a vertex, which voxel, and trilinear coordinates (alp) within the voxel.
    ind = floor((V - ll)./h)+[1,1,1];       % indices into BC(voxels)
    ind = ind(:,[2,1,3]); 
    % Take voxel coordinates and build indices into 8 voxel vertices (for each vertex)
    i = ind(:,1);
    j = ind(:,2);
    k = ind(:,3);
    DVix = [ sub2ind( sWp1, i  , j  , k  ), ...
             sub2ind( sWp1, i  , j  , k+1), ...
             sub2ind( sWp1, i  , j+1, k  ), ...
             sub2ind( sWp1, i  , j+1, k+1), ...
             sub2ind( sWp1, i+1, j  , k  ), ...
             sub2ind( sWp1, i+1, j  , k+1), ...
             sub2ind( sWp1, i+1, j+1, k  ), ...
             sub2ind( sWp1, i+1, j+1, k+1) ];
    SVix_hex = IM(DVix); % convert from voxel verties to hexahedral mesh vertices

    % Map each V to a tet(H_tet) and corresponding barycentric coordinates
    V2Htet = zeros(1, size(V,1));
    V2BaryCentric = zeros(size(V,1), 4);
    for ii=1:size(V,1)
        H_hex_idx = 0;
        for jj=1:size(H_hex,1)
            if all(SVix_hex(ii,:) == H_hex(jj,:))
                H_hex_idx = jj;
                break
            end
        end
        % H_hex_idx should be okay
        assert(H_hex_idx > 0);

        epsilon = 1e-8;
        % 6 candidate tets, H_cand_tets: 6x4
        H_cand_tets = H_tet(H_hex_idx*6-5:H_hex_idx*6,:);
        for jj=1:6
            % candidate tets vertices
            cand_TV = SV(H_cand_tets(jj,:),:);     % 4x3
            % try to compute the barycentric coordinate of V(ii,:) in cand_TV
            bcCoeff = barycentric_coordinates(V(ii,:),cand_TV(1,:),cand_TV(2,:),cand_TV(3,:),cand_TV(4,:));
            if all(bcCoeff + epsilon >= 0)
                V2Htet(ii) = H_hex_idx*6-6+jj;
                V2BaryCentric(ii,:) = bcCoeff;
                break
            end
        end

    end

    % build barycentric Jacobian matrix such that W * SV(:,k) = V(:,k)
    row = size(V,1);
    col = size(SV,1);
    W = sparse(row, col);
    for ii=1:size(V,1)
        W(ii,H_tet(V2Htet(ii),:)) = V2BaryCentric(ii,:);
    end

    V3(:,1) = W * SV(:,1);
    V3(:,2) = W * SV(:,2);
    V3(:,3) = W * SV(:,3);
    V = V3;
    T = H_tet;
%     plotTetrahedralMesh(SV,H_tet,'Style', '-b')
    
%     function [] = plotTetrahedralMesh( SV, H, varargin )
%             style = 'b:';
%             params_to_variables = containers.Map( {'Style'}, {'style'} );
%             v = 1;
%             while v <= numel(varargin)
%                 param_name = varargin{v};
%                 if isKey(params_to_variables,param_name)
%                   assert(v+1<=numel(varargin));
%                   v = v+1;
%                   % Trick: use feval on anonymous function to use assignin to this workspace
%                   feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
%                 else
%                   error('Unsupported parameter: %s',varargin{v});
%                 end
%                 v=v+1;
%             end
% 
%             bi = [ 1 2; 1 3; 1 4;  ...
%                    2 3; 2 4; 3 4];
%             lines = zeros( 1, 3*size(H,1)*size(bi,1) );
%             count = 1;
%             for a = 1:size(H,1)
%                 for n = 1:size(bi,1)
%                     lines( 1:3, count:count+2 ) = [ ...
%                         SV(H(a,bi(n,1)),1) SV(H(a,bi(n,2)),1), nan; ...
%                         SV(H(a,bi(n,1)),2) SV(H(a,bi(n,2)),2), nan; ...
%                         SV(H(a,bi(n,1)),3) SV(H(a,bi(n,2)),3), nan ];
%                     count = count + 3;
%                 end
%             end
%             plot3( lines(1,:), lines(2,:), lines(3,:), style);
%         end
end

