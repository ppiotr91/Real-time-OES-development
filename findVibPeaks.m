clc; clear; close all;
% For the first positive system, find the position of the band heads.
% Scatter plot these band heads on top of the spectrum generated in
% SpecAir.
%   ADDTIONAL NOTES ON USAGE:
%       -   Use the spectrums generated from specair with a very high
%           resolution. The data is called from .txt files and have the
%           following naming convenstion:
%               system_Tvibration_Trotation.txt
%           where, system = '0','1' or '2' corresponding to the first
%           negative, first positive and the second positive system
%           respectively. Tvibration and Trotation are the vibrational and
%           the rotational temperatures in Kelvin. 
%       
%       -   The wavelenghts for each transition are calcuated using
%           evalRotationalPeaks function. 
syst        = '1-';
V1          = 0;
V2          = 0;
J1          = 22;
J2          = 36;

[wlBH, ~]                   = evalBandHead(syst,V1,V2);
[wl, Bv1, J1, J2, txt]      = evalRotationalPeaks(syst,'P',V1,V2,J1,J2);
fortratDiagram(syst, V1, V2)

fName   = 'PWM-1000W-10-0.5Ar-10N2.csv';
[wavelength, intensity] = reshapeLF(fName);
intensity               = intensity/max(intensity);
[~,indexMax]            = max(intensity);
wavelength              = wavelength - (wavelength(indexMax) - wlBH);
inten                   = interp1(wavelength,intensity,wl);

figure
plot(wavelength, intensity);    hold on;
plot([wlBH, wlBH], [0,1],'--r');
scatter(wl(1:end-1),inten(1:end-1),'filled','DisplayName', 'First Negative System');
text(wl(1:end-1),inten(1:end-1),txt(1:end-1),'FontWeight','bold','FontSize',10,...
        'HorizontalAlignment','center');
axis([min(wl)-0.2 max(wl)+0.2 0 inf+0.3]);
xlabel('Wavelength (nm)');
ylabel('Normalized Intensity');

maxInten            = zeros(length(wl)-1,1);
maxWavelength       = zeros(length(wl)-1,1);
for i = 1:length(wl)-1
    indexLower                  = find(wavelength > wl(i+1), 1 );
    indexUpper                  = find(wavelength < wl(i), 1, 'last' );
    [~, indexMin]               = min(intensity(indexLower:indexUpper));
    indexMin                    = indexMin + indexLower - 1;
    [maxInten(i),indexMaxInten] = max(intensity(indexMin:indexUpper));
    indexMaxInten               = indexMaxInten + indexMin - 1;
    
    maxWavelength(i)    = wavelength(indexMaxInten);
end

scatter(maxWavelength,maxInten,'filled','DisplayName', 'First Negative System');

J1(end)     = [];
J2(end)     = [];
hc          = 1240 * 1e-7;      %eV - cm
g       = zeros(length(J1),1);
for i   = 1:length(J1)
    %even J2 case
    if (mod(J2(i),2) == 0)
        g(i)   = 2/3 * (J1(i) + J2(i) + 1);
         
    %odd J2 case
    else
        %       g(i)   = 1/3 * (J1(i) + J2(i) + 1);
        g(i) = 0;
        
    end
end

indexOddJ           = find(g == 0);
g(indexOddJ)        = [];
maxInten(indexOddJ) = [];
J1(indexOddJ)       = [];
J2(indexOddJ)       = [];

lhs     = log(maxInten ./ g);
rhs     = J1 .* (J1 + 1);
P       = polyfit(rhs, lhs, 1);
xfit    = linspace (min(rhs), max(rhs))';
yfit    = P(1) * xfit + P(2);

figure
scatter (rhs, lhs, 'filled', 'DisplayName', 'measured'); hold on;
plot (xfit, yfit, '--r', 'LineWidth', 2, 'DisplayName', ...
      sprintf('y = %.4E x + %.4G', P(1), P(2)));
xlabel('J"(J"+1)');
ylabel('ln(I/(J" + J"" + 1)');

legend('show');

tempRoteV   = -1 * Bv1 * hc /  P(1);
tempRotK    = tempRoteV / 8.617e-5;
txt2        = sprintf('T_r = %.3f eV = %.0f K', tempRoteV, tempRotK);
annotation('textbox', [.2 .3 0.9 0.0] , 'String', txt2 ,...
           'FitBoxToText', 'on');