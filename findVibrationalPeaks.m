clc; clear; close all;

% Plot the band heads of the various nitrogen systems (i.e. Find the
%   position of the rotationless vibrational trasistion). Plot as a function
%   of the up vibrational level vs wavelength position. 
%
%       PROCEDURE:
%           1. Define the upper and lower vibrational levels of interest.
%           2. Create a MeshGrid matrix to evaluate all possible
%           transtion between vibrational levels of the two electronic
%           systems.
%           3. Pass the meshGrid variables to roVib and obtain all the
%           wavelengths.
%               - NOTE:     roVib can allow for NxN vibrational or
%               rotational number inputs. However, the other will need to
%               be an integer number. 
%           4. Reshape the NxN matricies into a Mx1 matrix for scatter
%           plotting.
%           5. Create annotation lables for each transition.
%           6. Create a Scatter plot with labels.
%           7. Repeat process for any other system.

V1      = 0:4;
V2      = V1;
J1      = 0;
J2      = 0;
system  = '1-';

% IMPORT DATA
figure;
specData    = importdata('0_8000_8000.txt');
waveLength  = specData(:,1);
intensity   = specData(:,2)/max(specData(:,2));
plot(waveLength, intensity);    hold on;

% FIND BAND HEADS
% [V1,V2] = meshgrid(v1,v2);
% V1      = reshape(V1,[],1);
% V2      = reshape(V2,[],1);
txt     = cell(size(V1));
wlBH    = zeros(size(V1));

for i   = 1:length(V1)
    [wlBH(i), ~]    = evalBandHead( system, V1(i), V2(i) );
    txt(i)          = {sprintf('(%d, %d)',V1(i),V2(i))};
end

%SCATTER BAND HEADS
intenBH  = interp1(waveLength,intensity,wlBH);

scatter(wlBH,intenBH,'filled','DisplayName', 'First Negative System');
text(wlBH,intenBH+0.1,txt,'FontWeight','bold','FontSize',8,...
        'HorizontalAlignment','center');
axis([381 392 0 inf])


