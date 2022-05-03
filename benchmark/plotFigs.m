tdLittlePlanet = loadTDcellarray("littlePlanetTauERes_5001-25-2022_12-49.mat");
tdBeam = loadTDcellarray("beamTauERes_2001-23-2022_21-08.mat");
tdCantileverPull = loadTDcellarray("tauEResPull_5001-25-2022_12-27.mat");
% tdCantileverPull = loadTDcellarray("tauEResPull_5001-24-2022_16-41.mat");

ERresLittlePlanet = AnalyseTDtauE(tdLittlePlanet,logspace(-10,0,50),2);
ERresBeam = AnalyseTDtauE(tdBeam,logspace(-10,0,20),3);
ERresCantileverPull = AnalyseTDtauE(tdCantileverPull,logspace(-10,0,50),2);

allErrors = [ERresLittlePlanet(:,2) ; ERresBeam(:,2) ; ERresCantileverPull(:,2)];
% allErrors = [ERresCantileverPull(:,2)];
leftYlims = minmax(log10(allErrors'));

[f1,f2] = plotERres(ERresLittlePlanet, leftYlims);
saveas(f1, "benchmark/autogenFigs/littleplanetPareto.pdf");
saveas(f2, "benchmark/autogenFigs/littleplanet.pdf");
% 
[f1,f2] = plotERres(ERresBeam, leftYlims);
saveas(f1, "benchmark/autogenFigs/beamPareto.pdf");
saveas(f2, "benchmark/autogenFigs/beam.pdf");
% 
[f1,f2] = plotERres(ERresCantileverPull, leftYlims);
saveas(f1, "benchmark/autogenFigs/cantileverPullPareto.pdf");
saveas(f2, "benchmark/autogenFigs/cantileverPull.pdf");

function  [f1, f2]=plotERres(ERres, leftYLims)
    speedup = ERres(:,1);
    maxError = ERres(:,2);
    tre = ERres(:,3);
    logtre = log10(tre);

    figure;
    scatter(1./speedup,log10(maxError));
    title('Pareto Speedup');
    ylabel('log maxError');
    xlabel('speed factor');
    f1 = gca;
    
    figure;
    yyaxis left;
    plot(logtre, log10(maxError));
    ylabel('log max error (meters) size scaled');
%     ylim(leftYLims);
    yyaxis right;

    plot(logtre,1./speedup);
    title('Impact of the elastic threshold');
    xlabel('log tau E');
    ylabel('factor');
    f2 = gca;
end
