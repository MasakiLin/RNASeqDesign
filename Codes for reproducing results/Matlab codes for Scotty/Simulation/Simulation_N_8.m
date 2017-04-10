%repress figure
figure('visible','on');

set(gcf,'Visible','off')              % turns current figure "off"
set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"

files = dir('Scotty.dispersion.5.N.8.*.txt');
myArray = zeros(149,10,100)
for i =1:20
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '8', '8', '1', '1.4', '0.0001' ,'50', '0', '0', '0', 'Inf', '150', '1000000', '20000000', '100', '100', '100', 'result');
end

i=1
dlmwrite('Simulation_dispersion_5_Scotty_N_8.txt',myArray(:,:,i))
for i = 2:20
dlmwrite('Simulation_dispersion_5_Scotty_N_8.txt',myArray(:,:,i), '-append')
end

%scottyEstimate( fileName, nControlSamples, nTestSamples, outputTag, ...
%    fc, pCut, minPercDetected, costPerRepControl, costPerRepTest, costPerMillionReads, totalBudget, ...
%    maxReps, minReadsPerRep, maxReadsPerRep, minPercUnBiasedGenes, pwrBiasCutoff, alignmentRate, ...
%    outputDirectory)
