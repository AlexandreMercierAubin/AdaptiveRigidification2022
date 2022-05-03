function plotBenchmarkCSV(fileName)
    f = figure();
    M = readtable(fileName);
    
    adaptiveLineX = [];
    adaptiveLineY = [];
    lineX = [];
    lineY = [];
    for i = 1:size(M,1)
        name = strrep(cell2mat(M{i,1}),'beam','');
        if strcmp(cell2mat(M{i,2}), 'true')
            adaptiveLineY = [adaptiveLineY, M{i,3}];
            adaptiveLineX = [adaptiveLineX, str2double(name)];
        else
            lineY = [lineY, M{i,3}];
            lineX = [lineX, str2double(name)];
        end
    end
    hold on;
    title("adaptive rigidification vs normal")
    semilogy(adaptiveLineX, adaptiveLineY, 'Color', [1,0.25,0.25]);
    semilogy(lineX, lineY, 'Color', [0,0,1]);
end

