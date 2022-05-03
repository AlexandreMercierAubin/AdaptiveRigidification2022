function customCMap = makeCmap(data, middle,topColor, indexColor, bottomColor)
        %https://www.mathworks.com/matlabcentral/answers/305073-colormap-fixed-middle-value
        L = numel(data);
        
        % Calculate where proportionally indexValue lies between minimum and
        % maximum values
        largest = max(max(data));
        smallest = min(min(data));
        index = L*abs(middle-smallest)/(largest-smallest);
        % Create color map ranging from bottom color to index color
        % Multipling number of points by 100 adds more resolution
        customCMap1 = [linspace(bottomColor(1),indexColor(1),100*index)',...
                    linspace(bottomColor(2),indexColor(2),100*index)',...
                    linspace(bottomColor(3),indexColor(3),100*index)'];
        % Create color map ranging from index color to top color
        % Multipling number of points by 100 adds more resolution
        customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
                    linspace(indexColor(2),topColor(2),100*(L-index))',...
                    linspace(indexColor(3),topColor(3),100*(L-index))'];
        customCMap = [customCMap1;customCMap2];  % Combine colormaps
end