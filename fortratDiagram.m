function [] = fortratDiagram(system, V1, V2, numJ)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5
    numJ    = 50;
end

% Evaluate P-branch
J1i     = 0;
J1f     = numJ;

[wl, ~, J2, ~] = evalRotationalPeaks(system, 'P', V1, V2, J1i, J1f);
scatter (wl, J2,'r','filled','DisplayName','P-Branch'); hold on;

% Evaluate R-branch
J1i     = 1;
J1f     = 1 + numJ;

[wl, ~, J2, ~] = evalRotationalPeaks(system, 'R', V1, V2, J1i, J1f);
scatter (wl, J2,'k','filled','DisplayName','R-Branch');

xlabel('Wavelength (nm)');
ylabel('Rotational quantum number, J''''');
title(sprintf('Fortrat Diagram for (v'', v'''') = (%d,%d)',V1, V2));
legend('show')

end

