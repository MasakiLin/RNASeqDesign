%repress figure
figure('visible','on');

set(gcf,'Visible','off')              % turns current figure "off"
set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"

files = dir('Data.2*.txt');
myArray = zeros(149,10,100)
for i =1:10
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '2', '2', '1', '1.2', '0.05' ,'50', '0', '0', '0', 'Inf', '150', '1000000', '6700000', '100', '100', '100', 'result');
end

i=1
dlmwrite('Result_n_2.txt',myArray(:,:,i))
for i = 2:10
dlmwrite('Result_n_2.txt',myArray(:,:,i), '-append')
end 

files = dir('Data.4.*.txt');
myArray = zeros(149,10,100)
for i =1:10
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '4', '4', '1', '1.2', '0.05' ,'50', '0', '0', '0', 'Inf', '150', '1000000', '6700000', '100', '100', '100', 'result');
end

i=1
dlmwrite('Result_n_4.txt',myArray(:,:,i))
for i = 2:10
dlmwrite('Result_n_4.txt',myArray(:,:,i), '-append')
end 


files = dir('Data.6.*.txt');
myArray = zeros(149,10,100)
for i =1:10
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '6', '6', '1', '1.2', '0.05' ,'50', '0', '0', '0', 'Inf', '150', '1000000', '6700000', '100', '100', '100', 'result');
end

i=1
dlmwrite('Result_n_6.txt',myArray(:,:,i))
for i = 2:10
dlmwrite('Result_n_6.txt',myArray(:,:,i), '-append')
end 

files = dir('Data.8.*.txt');
myArray = zeros(149,10,100)
for i =1:10
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '8', '8', '1', '1.2', '0.05' ,'50', '0', '0', '0', 'Inf', '150', '1000000', '6700000', '100', '100', '100', 'result');
end


i=1
dlmwrite('Result_n_8.txt',myArray(:,:,i))
for i = 2:10
dlmwrite('Result_n_8.txt',myArray(:,:,i), '-append')
end 



files = dir('Data.10.*.txt');
myArray = zeros(149,10,100)
for i =1:10
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '10', '10', '1', '1.2', '0.05' ,'50', '0', '0', '0', 'Inf', '150', '1000000', '6700000', '100', '100', '100', 'result');
end

i=1
dlmwrite('Result_n_10.txt',myArray(:,:,i))
for i = 2:10
dlmwrite('Result_n_10.txt',myArray(:,:,i), '-append')
end 


files = dir('Data.12.*.txt');
myArray = zeros(149,10,100)
tmp=files(1).name
myArray(:,:,1)=scottyEstimate(tmp, '12', '12', '1', '1.2', '0.05' ,'50', '0', '0', '0', 'Inf', '150', '1000000', '6700000', '100', '100', '100', 'result');

i=1
dlmwrite('Result_n_12.txt',myArray(:,:,i))


%%%%%%%%%%%%%%%%%%%%%%% Depth

%repress figure
figure('visible','on');

set(gcf,'Visible','off')              % turns current figure "off"
set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"

files = dir('Data.2*.txt');

Data=textread('p.cut.Data.2.txt','%f')
myArray = zeros(149,10,100)
for i =1:50
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '2', '2', '1', '1.2', num2str(Data(i)) ,'50', '0', '0', '0', 'Inf', '150', '1000000', '8000000', '100', '100', '100', 'result');
end

i=1
dlmwrite('Result_n_2_Depth.txt',myArray(:,:,i))
for i = 2:50
dlmwrite('Result_n_2_Depth.txt',myArray(:,:,i), '-append')
end 

files = dir('Data.4.*.txt');
Data=textread('p.cut.Data.4.txt','%f')
myArray = zeros(149,10,100)
for i =1:50
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '4', '4', '1', '1.2', num2str(Data(i)) ,'50', '0', '0', '0', 'Inf', '150', '1000000', '8000000', '100', '100', '100', 'result');
end

i=1
dlmwrite('Result_n_4_Depth.txt',myArray(:,:,i))
for i = 2:50
dlmwrite('Result_n_4_Depth.txt',myArray(:,:,i), '-append')
end 

files = dir('Data.8.*.txt');
Data=textread('p.cut.Data.8.txt','%f')
myArray = zeros(149,10,100)
for i =29:50
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '8', '8', '1', '1.2', num2str(Data(i)) ,'50', '0', '0', '0', 'Inf', '150', '1000000', '8000000', '100', '100', '100', 'result');
end

i=1
dlmwrite('Result_n_8_Depth.txt',myArray(:,:,i))
for i = 2:50
dlmwrite('Result_n_8_Depth.txt',myArray(:,:,i), '-append')
end 


files = dir('Data.16.*.txt');
Data=textread('p.cut.Data.16.txt','%f')
myArray = zeros(149,10,100)
for i =1:50
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '16', '16', '1', '1.2', num2str(Data(i)) ,'50', '0', '0', '0', 'Inf', '150', '1000000', '8000000', '100', '100', '100', 'result');
end

i=1
dlmwrite('Result_n_16_Depth.txt',myArray(:,:,i))
for i = 2:50
dlmwrite('Result_n_16_Depth.txt',myArray(:,:,i), '-append')
end 
