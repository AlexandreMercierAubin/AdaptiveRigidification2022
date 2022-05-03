%takes in 
function readTDcsv(fileNames, maxFrames)
    if nargin < 2
        maxFrames = -1;
    end
    fontSize = 8;
    fontName = 'Linux Biolinum O';
    hfig = figure('Renderer', 'painters', 'Units','Inches', 'Position', [0,0,3.3125,2]);%
    set(gcf,'color','w');
    hold on;
    %for i = 1:numel(fileNames)
       myTable1 = readtable(fileNames(1));
       myTable2 = readtable(fileNames(2));
       if maxFrames > 0
           myTable1 = myTable1(1:maxFrames, :);
           myTable2 = myTable2(1:maxFrames, :);
       end
       %col 10 is the total sim time and 11 is the time to compute the
       %energy
       times1 = table2array(myTable1(:,10)) - table2array(myTable1(:,11));
       st1 = sum(times1)
       times2 = table2array(myTable2(:,10)) - table2array(myTable2(:,11));
       st2 = sum(times2)
       timeNoContact1 = times1 - table2array(myTable1(:,15));
       stn1 = sum(timeNoContact1)
       timeNoContact2 = times2 - table2array(myTable2(:,15));
       stn2 = sum(timeNoContact2)

       totSpeedup = st2/st1
       totSpeedupNoContact = stn2/stn1

       speedup = times2./times1;
       meanspeed = mean(speedup)

       speedupNoContact = timeNoContact2./timeNoContact1;
       meanspeedNoContact = mean(speedupNoContact)

       plot([1:size(times1,1)]', speedup);
       plot([1, size(times1,1)], [meanspeed, meanspeed], 'g');
       string(num2str(speedup));
       ylabel('Speedup Factor');
       xlabel('Time Step');
    %end
    
    handle = gca;
    handle.Color = [0.95, 0.95, 0.95];
    handle.XLim = [0, size(times1,1)];
    handle.TickLength = [0,0];
    handle.YGrid = 'on';
    handle.GridColor = [1 1 1];
    handle.XColor = [0 0 0];
    set(gca, 'fontSize', fontSize, 'fontName', fontName);
%     set(hfig,'Units','Inches');
    pos = get(hfig,'Position');
    set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperPosition',[0,0,3.3125,2]);%
    printSvg = strrep(fileNames(1),'.csv','.svg');
%     exportgraphics(gca,"benchmark/"+printPdf,'ContentType','vector');
    saveas(gca,"benchmark/"+printSvg)
    %legend(fileNames);
    
    
end

