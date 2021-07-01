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
%           evalRotaitonalPeaks function. 

specData    = importdata('0_2000_8000.txt');
waveLength  = specData(:,1) + 0.1184;
intensity   = specData(:,2)/max(specData(:,2));
plot(waveLength, intensity);    hold on;

[wl, J1,J2, txt]  = evalRotationalPeaks('1-','P',0,0,0,40);
wl  = wl + 0.004;
inten   = interp1(waveLength,intensity,wl);

scatter(wl,inten,'filled','DisplayName', 'First Negative System');
text(wl,inten+0.1,txt,'FontWeight','bold','FontSize',8,...
        'HorizontalAlignment','center');
    
axis([300 550 0 inf]);



% v1  = 0:3;
% v2  = 0:6;
% J1  = 0;
% J2  = 0;
% [V1,V2] = meshgrid(v1,v2);
% wl      = roVib('1-',V1, V2, J1, J2, 0);
% 
% V1      = reshape(V1,[],1);
% V2      = reshape(V2,[],1);
% wl      = reshape(wl,[],1);

% 
% txt     = cell(size(V1));
% for i   = 1:length(V1)
%     txt(i)  = {sprintf('(%d, %d)',V1(i),V2(i))};
% end
