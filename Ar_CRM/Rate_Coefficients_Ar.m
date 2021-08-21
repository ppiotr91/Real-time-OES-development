clear
clc
close all

[~,sheets] = xlsfinfo('Q_e_Ar.xlsx');

data = [];
data_curr = [];
for i = 1:numel(sheets)
    if i == 1
        data = readtable('Q_e_Ar.xlsx', 'Sheet' ,i);
    else
        data_curr = readtable('Q_e_Ar.xlsx', 'Sheet' ,i);
    end
    data = [data; data_curr];
end

Te = linspace(0.1, 10^3, 10^4); %eV
Tg = 5000; %Kelvin

