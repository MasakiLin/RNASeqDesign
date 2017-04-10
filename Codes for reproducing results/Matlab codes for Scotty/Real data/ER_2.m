%repress figure
figure('visible','on');

set(gcf,'Visible','off')              % turns current figure "off"
set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"

files = dir('ER.Data.2.*.txt');
myArray = zeros(149,10,100)
for i =1:10
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '2', '2', '1', '1.4', '0.0001' ,'50', '0', '0', '0', 'Inf', '150', '1000000', '95000000', '100', '100', '100', 'result');
end

i=1
dlmwrite('ER_result_n_2.txt',myArray(:,:,i))
for i = 2:10
dlmwrite('ER_result_n_2.txt',myArray(:,:,i), '-append')
end

%scottyEstimate( fileName, nControlSamples, nTestSamples, outputTag, ...
%    fc, pCut, minPercDetected, costPerRepControl, costPerRepTest, costPerMillionReads, totalBudget, ...
%    maxReps, minReadsPerRep, maxReadsPerRep, minPercUnBiasedGenes, pwrBiasCutoff, alignmentRate, ...
%    outputDirectory)
