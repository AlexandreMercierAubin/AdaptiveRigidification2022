function ERres = AnalyseTDtauE(td, x, dim)
    %TD: timing datas
    %number of test entries in total including the elastic sim last
    %x list rigid thresholds
    %y list elastic multiplyers
    %elastic positions in time
    
    testi = numel(td);
    for frame = 1:size(td{testi}.p,2) 
        pb{frame} = reshape(td{testi}.p(:,frame), dim, [])';
    end
    sceneLength = max(pdist(pb{frame}));

    maxErr = zeros(testi-1,1);
    speedup = zeros(testi-1,1);
    timeEl = sum(td{testi}.log(2,:));
    tre = zeros(testi-1,1);
    trr = zeros(testi-1,1);
    sumDofs = zeros(testi-1,1);
    harmonicMeanSpeedup = zeros(testi-1,1);

    for i = 1:testi-1
        for frame = 1:size(td{i}.p,2)
            pi = reshape(td{i}.p(:,frame), dim, [])';
    %         errDistance = sum(sqrt(sum((pb{frame}-pi).^2,2)));
            errDistance = max(sqrt(sum((pb{frame}-pi).^2,2)));
            maxErrTmp = max(errDistance) / sceneLength;
            maxErr(i) = max(maxErr(i),maxErrTmp);
        end
        time = sum(td{i}.log(2,:));
        speedup(i) = timeEl/time;
        sumDofs(i) = sum(td{i}.logCounts(7,:)); %7 is the id for totaldofs in td
        perFrameSpeedup = td{i}.log(2,:)./td{testi}.log(2,:);
        harmonicMeanSpeedup(i) = harmmean(perFrameSpeedup);
        
        tre(i) = x(i);
        trr(i) = tre(i)*0.1;
    end
    dofsRatio = sumDofs./sum(td{testi}.logCounts(7,:));
    ERres = [speedup,maxErr,tre,trr,sumDofs,dofsRatio];

end


