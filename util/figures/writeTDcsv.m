function writeTDcsv(td, testname, names)
% If you want to be very thorough, you can even install libertine and generate plots that have fonts that match the paper.
% 
% 1) Download libertine from sourceforge
% 2) Right click and install the fonts 
%     - Copying to windows/fonts folder doesn't seem to be enough
%     - Might only need LinLibertine_Rah
% 3) Use fontName = 'Linux Libertine' in matlab
% 4) Use print command with -dsvg option
% 5) Open in inkscape
% 6) Save as pdf with default options
% 
% The inkscape thing....  
% 
% n Inkscape:
% - open, clip to crop box (other default settings are probably fine)
% - save as, (pdf, default setting fine, i.e., embed fonts)
% 
% that will get the fonts embedded.
% 
% Otherwise... here are the basics for setting up a matlab script for generating a nice plot:
% 
% fontSize = 10;
% fontName = 'Times New Roman';
% hfig = figure('Renderer', 'painters', 'Position', [10 10 600 400])
% set(gcf,'color','w');
% hold on;
% 
% % do your plotting 
% 
% set(gca, 'fontsize', fontSize, 'fontname', fontName);
% set(hfig,'Units','Inches');
% pos = get(hfig,'Position');
% set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(hfig,'figures/figRMagAllLambdaLinear','-dpdf','-r0')
% 
% Adjust the parameters as necessary... likewise, when doing your plotting, keep in mind consistency, and how things show up in the doc (i.e., linewidth, colour choices, etc)
%WRITETDCSV Summary of this function goes here
%   Detailed explanation goes here
    assert(numel(td) == numel(names));

    for i = 1:numel(td)
        testFile = testname + names(i) + ".csv";
        hc = td{i}.labelsCountData;
        h = td{i}.labels;
        headers = sprintf('%s,', hc);
        headers2 = sprintf('%s,', h);
        headers2 = headers2(1:end-1);
        headers =[headers, headers2, '\n'];

        fileID = fopen(testFile,'w');
        nbytes = fprintf(fileID,headers,1);
        for j = 1:size(td{i}.log,2);
            row = [td{i}.logCounts(:,j);td{i}.log(:,j)];
            rowS = sprintf('%d,', row);
            rowS = [rowS(1:end-1), '\n'];
            nbytes = fprintf(fileID,rowS,1);
        end
        fclose(fileID);
    end
end

