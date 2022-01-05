function [ WJH, J2BH ] = evalBandHead( system, V1, V2 )
% Find the band head position for a given transition between V' and V".
% Bandhead wavelenght is taken to be the wavelenght as the maximum wavelenght 
% on a P-Branch parabola. 
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
%
%       OUTPUTS:-
%           - WJH:      Wavelenght position of the band head (in nm)
%           - J2BH:     Rotational quantum number of the lower level. 

[wl, ~, J2, ~]      = evalRotationalPeaks(system, 'P', V1, V2, 0, 50);
[WJH, i]            = max(wl);
J2BH                = J2(i);

end

