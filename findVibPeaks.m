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
J1          = 20;
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
scatter(wl,inten,'filled','DisplayName', 'First Negative System');
text(wl,inten,txt,'FontWeight','bold','FontSize',10,...
        'HorizontalAlignment','center');
axis([min(wl)-0.2 max(wl)+0.2 0 inf+0.3]);

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

% [maxInten,indexMaxInten]    = findpeaks(intensity(indexLower:indexUpper));
% maxWavelength(i)            = wavelength(indexLower + indexMaxInten - 1);





