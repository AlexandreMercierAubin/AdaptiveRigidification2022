% M taken from slack (loaded from file?)
function plotTauData(data, M, N) 

    speedup  = reshape( data(:,1), M, N )';
    error = reshape( data(:,2), M, N )';
    relErr = error ./ max(max(error));
    dofs = reshape( data(:,5), M, N )';
    dofsRatio = reshape( data(:,6), M, N )';
    logError = reshape( log10( data(:,2) ), M, N )';
    logRatio = reshape( log10( data(:,3) ), M, N )';
    logtauR  = reshape( log10( data(:,4) ), M, N )';

    invErr = 1./error;
    lie = log10(invErr);
    lre = log10(relErr);
    
    factor = 1./speedup;
    logFactor = log10( factor );
    cMap = interp1([0;1],[1 0 0; 0 1 0],linspace(0,1,256));
    cMapInv = interp1([0;1],[0 0 0; 1 0 0],linspace(0,1,256));
    
    clf;
    i = 0;
    %i=i+1;doSubplot( i, relErr, 'relError' );
    i=i+1;doSubplot( i, -lre, '-log relError' , M, N,cMap);
    customCmap = makeCmap(speedup-1,0, [0,1,0,0.2], [1,1,1,0.2], [1,0,0,0.2]);
    i=i+1;doSubplot( i, speedup-1, 'speedup-1' , M, N,customCmap);
    i=i+1;doSubplot( i, factor, 'factor (inv speedup)' , M, N,cMap);
    %i=i+1;doSubplot( i, logFactor, 'log Factor' );
    %i=i+1;doSubplot( i, logError, 'log error' );
    i=i+1;doSubplot( i, -(factor + lre/10), 'score4' , M, N,cMap);
    i=i+1;doSubplot( i, 1./(1./logFactor + 3./lre), '1/(1/log inv speedup + 3/relError)', M, N,cMap );
    i=i+1;doSubplot( i, 1./(1./factor + 1./lre), '1/(speedup+1/relError)' , M, N,cMap);
    i=i+1;doSubplot( i, -lre .* (speedup-1), 'log relError(speedup-1)' , M, N,cMap);
    %     i=i+1;doSubplot( i, dofs, 'dofs' , M, N,cMapInv);
    
    i=i+1;doSubplot( i, log10(error), 'relErr' , M, N,cMapInv);
    hold on;
    alpha(0.5);
    doSubplot( i, speedup-1, 'speedup-1 & relErr' , M, N,customCmap);
    alpha(0.5);
    
%     i=i+1;doSubplot( i, harmonic, 'harmonic' , M, N,cMap);

    customCmap = makeCmap(dofsRatio,1, [1,0,0,0.2], [1,1,1,0.2], [0,1,0,0.2]);
    i=i+1;doSubplot( i, dofsRatio, 'dofs ratio' , M, N,customCmap);
    
    function doSubplot( sp, val, name, M, N, cMap)
        x = linspace(-10,0, M);
        y = linspace(-3,0, N);
        [X,Y] = meshgrid(x,y);
        ax1 = subplot(3,3,sp);
        fh = imagesc( x, y, val ); 
        axis xy;
        colormap(gca,cMap);
        colorbar;
        hold on;
        numContours = 7;
        [~,c] = contour( X, Y, val, numContours, 'k' ); 
        c.LineWidth = 1;
        axis equal
        axis([min(x),max(x),min(y),max(y)]);
        title(name);
        xlabel('log tau E');
        ylabel('ratio');
    end

    
end


    
