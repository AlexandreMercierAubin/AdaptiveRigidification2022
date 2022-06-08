function [AdjagencyMatrix,TetsPerParticle,Graph] = elementAdjacencyMatrix(T,N)
    TetsPerParticle = cell(N, 1, 1);
    Tsorted = sort(T,2);
    sharedFaces = containers.Map('KeyType','char','ValueType','any');
    iimat = [];
    jjmat = [];
    vals = [];
    for i = 1:1:size(Tsorted, 1)
        for j = 1:size(Tsorted,2)
                TetsPerParticle{Tsorted(i,j)} = [TetsPerParticle{Tsorted(i,j)},i];
        end

        if size(Tsorted,2) == 3
            % Check all triangles with higher index, to see which
            % indices are shared with our triangle... and if exactly
            % two are shared then we share an edge with that triangle!
            % That is, we've found the adjacent triangles.
            edge1 = key_enc([Tsorted(i,1),Tsorted(i,2)]);
            edge2 = key_enc([Tsorted(i,1),Tsorted(i,3)]);
            edge3 = key_enc([Tsorted(i,2),Tsorted(i,3)]);
            edges = [edge1,edge2,edge3];
            for j = 1:numel(edges)
                if isKey(sharedFaces,edges{j})
                    matching = sharedFaces(edges{j});
                    iimat = [iimat;i;matching];
                    jjmat = [jjmat;matching;i];
                    vals = [vals;1;1];
                else
                    sharedFaces(edges{j}) = i;
                end
            end
        
        else
%             TetsPerParticle(inds) = cellfun(@(x) [x, i], TetsPerParticle(inds), 'UniformOutput', false);

            face1 = key_enc([Tsorted(i,1),Tsorted(i,2),Tsorted(i,3)]);
            face2 = key_enc([Tsorted(i,1),Tsorted(i,2),Tsorted(i,4)]);
            face3 = key_enc([Tsorted(i,1),Tsorted(i,3),Tsorted(i,4)]);
            face4 = key_enc([Tsorted(i,2),Tsorted(i,3),Tsorted(i,4)]);
            faces = [face1,face2,face3,face4];
            for j = 1:numel(faces)
                if isKey(sharedFaces,faces{j})
                    matching = sharedFaces(faces{j});
                    iimat = [iimat;i;matching];
                    jjmat = [jjmat;matching;i];
                    vals = [vals;1;1];
                else
                    sharedFaces(faces{j}) = i;
                end
            end
        end
    end
    AdjagencyMatrix = sparse(iimat,jjmat,vals,size(T,1),size(T,1));
    Graph = graph(AdjagencyMatrix);

    function key1_char =  key_enc(key1)
        [m,n] = size(key1);
        key1_char=repmat({char(ones(1,16*n))},m,1);
        tmp = num2hex(key1);
        for ii=1:m  
        key1_char{ii,:} = reshape((tmp(ii:m:end,:).'),1,16*n);
        end
    end
end

