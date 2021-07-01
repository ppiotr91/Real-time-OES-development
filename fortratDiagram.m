function [] = fortratDiagram(system, V1, V2, numJ)
% Generates fortrat diagrams for the various nitrogen systems.
%   ADDITIONAL NOTES ON USAGE:
%       INPUTS:- 
%           - system:   System associated with transition between
%                           electronic states. Choose either of the
%                           following:
%                               1. '1+' First positive system.
%                               2. '2+' Second positive system.
%                               3. '1-' First negative system.
%           - V1:       Vibrational level of emitting state (integer array)
%           - V2:       Vibrational level of lower state (integer array)
%           - numJ(OPTIONAL):
%                       Number of transitions of each branch to plot.
%                       This argument is optional. Default value is 50.
%
%       OUTPUTS:-
%           - No outputs but generates a figure.
%

if nargin < 5
    numJ    = 50;
end

% Evaluate P-branch
J1i     = 0;
J1f     = numJ;

figure
[wl, ~, J2, ~] = evalRotationalPeaks(system, 'P', V1, V2, J1i, J1f);
scatter (wl, J2,'r','filled','DisplayName','P-Branch'); hold on;

[wlBandhead, J2Bandhead]    = evalBandHead(system,V1,V2);
txt                         = sprintf('(%d, %d)\n %.2f nm',J2Bandhead-1, J2Bandhead, wlBandhead);
scatter (wlBandhead, J2Bandhead, 'b', 'filled', 'DisplayName', 'Bandhead');

% Evaluate R-branch
J1i     = 1;
J1f     = 1 + numJ;

[wl, ~, J2, ~] = evalRotationalPeaks(system, 'R', V1, V2, J1i, J1f);
scatter (wl, J2,'k','filled','DisplayName','R-Branch');

xlabel('Wavelength (nm)');
ylabel('Rotational quantum number, J''''');
title(sprintf('Fortrat Diagram for (v'', v'''') = (%d,%d)',V1, V2));
legend('show')
axis([-inf wlBandhead+2.2 -inf inf]);
text(wlBandhead+1.2,J2Bandhead,txt,'FontWeight','bold','FontSize',12,...
        'HorizontalAlignment','center');

end

