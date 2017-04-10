%repress figure
figure('visible','on');

set(gcf,'Visible','off')              % turns current figure "off"
set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"

files = dir('Data.5*.txt');
myArray = zeros(149,10,100)
for i =1:100
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '5', '5', '1', '1.15', '0.01', '50', '0', '0', '0', 'Inf', '150', '200000', '2000000', '100', '100', '100', 'result');
end

i=1
dlmwrite('Result_n_5.txt',myArray(:,:,i))
for i = 2:100
dlmwrite('Result_n_5.txt',myArray(:,:,i), '-append')
end 

files = dir('Data.10.*.txt');
myArray = zeros(149,10,100)
for i =1:100
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '10', '10', '1', '1.15', '0.01', '50', '0', '0', '0', 'Inf', '150', '200000', '2000000', '100', '100', '100', 'result');
end

i=1
dlmwrite('Result_n_10.txt',myArray(:,:,i))
for i = 2:100
dlmwrite('Result_n_10.txt',myArray(:,:,i), '-append')
end 

files = dir('Data.20.*.txt');
myArray = zeros(149,10,100)
for i =1:100
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '20', '20', '1', '1.15', '0.01', '50', '0', '0', '0', 'Inf', '150', '200000', '2000000', '100', '100', '100', 'result');
end

i=1
dlmwrite('Result_n_20.txt',myArray(:,:,i))
for i = 2:100
dlmwrite('Result_n_20.txt',myArray(:,:,i), '-append')
end 


files = dir('Data.30.*.txt');
myArray = zeros(149,10,100)
for i =1:100
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '30', '30', '1', '1.15', '0.01', '50', '0', '0', '0', 'Inf', '150', '200000', '2000000', '100', '100', '100', 'result');
end

i=1
dlmwrite('Result_n_30.txt',myArray(:,:,i))
for i = 2:100
dlmwrite('Result_n_30.txt',myArray(:,:,i), '-append')
end 




%%%%%%%%%%%%%%%%%%
%repress figure
figure('visible','on');

set(gcf,'Visible','off')              % turns current figure "off"
set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"

files = dir('Data.5.stronger*.txt');
myArray = zeros(149,10,100)
for i =1:100
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '5', '5', '1', '2', '0.05', '50', '0', '0', '0', 'Inf', '150', '200000', '2000000', '100', '100', '100', 'result');
end

i=1
dlmwrite('Result_n_5_stronger_new.txt',myArray(:,:,i))
for i = 2:100
dlmwrite('Result_n_5_stronger_new.txt',myArray(:,:,i), '-append')
end 


figure('visible','on');

set(gcf,'Visible','off')              % turns current figure "off"
set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"

files = dir('Data.20*.txt');
myArray = zeros(149,10,100)
for i =1:100
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '20', '20', '1', '2', '0.05', '50', '0', '0', '0', 'Inf', '150', '200000', '2000000', '100', '100', '100', 'result');
end

i=1
dlmwrite('Result_n_20_stronger_new.txt',myArray(:,:,i))
for i = 2:100
dlmwrite('Result_n_20_stronger_new.txt',myArray(:,:,i), '-append')
end 



figure('visible','on');

set(gcf,'Visible','off')              % turns current figure "off"
set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"

files = dir('Data.30*.txt');
myArray = zeros(149,10,100)
for i =1:100
    tmp=files(i).name
    myArray(:,:,i)=scottyEstimate(tmp, '30', '30', '1', '2', '0.05', '50', '0', '0', '0', 'Inf', '150', '200000', '2000000', '100', '100', '100', 'result');
end

i=1
dlmwrite('Result_n_30_stronger_new.txt',myArray(:,:,i))
for i = 2:100
dlmwrite('Result_n_30_stronger_new.txt',myArray(:,:,i), '-append')
end 

