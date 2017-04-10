%repress figure
figure('visible','on');

set(gcf,'Visible','off')              % turns current figure "off"
set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"

files = dir('Stage.Data.6.*.txt');
myArray = zeros(149,10,100)
for i =1:10
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '6', '6', '1', '1.4', '0.0001' ,'50', '0', '0', '0', 'Inf', '150', '1000000', '95000000', '100', '100', '100', 'result');
end

i=1
dlmwrite('Stage_result_n_6.txt',myArray(:,:,i))
for i = 2:10
dlmwrite('Stage_result_n_6.txt',myArray(:,:,i), '-append')
end 
%scottyEstimate( fileName, nControlSamples, nTestSamples, outputTag, ...
%    fc, pCut, minPercDetected, costPerRepControl, costPerRepTest, costPerMillionReads, totalBudget, ...
%    maxReps, minReadsPerRep, maxReadsPerRep, minPercUnBiasedGenes, pwrBiasCutoff, alignmentRate, ...
%    outputDirectory)