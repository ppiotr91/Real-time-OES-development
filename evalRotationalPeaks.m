function [WAVELENGTH,J1,J2,TXT] = evalRotationalPeaks(system, branch, V1, V2, J1i, J1f)
% Use this function to evaluate the wavelength of the rotational lines
% between upper level rotational quantum numbers, J1, for any vibrational
% system. The P-Branch will only provide outputs with a delJ of +1.
% R-Branch will provide values for delJ of -1.
%   ADDITIONAL NOTES ON USAGE:
%       INPUTS:-
%           - system:   System associated with transition between
%                           electronic states. Choose either of the
%                           following:
%                               1. '1+' First positive system.
%                               2. '2+' Second positive system.
%                               3. '1-' First negative system.
%           - branch:   Branch associated with the change in the
%                           rotataional quantum number J. Choose either of
%                           the following.
%                               1. 'P' P-Branch (delJ = +1)
%                               2. 'R' R-Branch (delJ = -1)
%           - V1:       Vibrational level of emitting state (integer array)
%           - V2:       Vibrational level of lower state (integer array)
%           - J1i:      Lower bound on the upper rotational transition.
%               (integer)
%           - J1f:      Upper bound on the upper rotational transition.
%               (integer)
%
%       OUTPUTS:-
%           - WAVELENGTH: Wavelength of the transitions of interest in nm.
%           - J1: Matrix containing all the upper level rotational quantum
%               numbers.
%           - J2: Matrix containing all the lower level rotational quantum
%               numbers.
%           - TXT: A cell array with a summary of each transition in form
%           of '(J1,J2). Can be used for plot annotation. 

switch branch
    %For P-Branch, delJ = +1. J1 can be evaluated upto 0.
    case 'P'
        J1  = J1i:J1f;
        J2  = J1 + 1;
        
    %For R-Branch, delJ = -1. J1 CANNOT be evaluated at 0.     
    case 'R'
        if J1i == 0
            warning('Upper Rotational quantum number cannot be evaluated at 0');
            J1i     = 1;
        end
        J1  = J1i:J1f;
        J2  = J1 - 1;
    otherwise
        error('Invalid input for branch selected. Choose "R" or "P" branches');
end

WAVELENGTH = roVib(system,V1,V2,J1,J2,0);

TXT     = cell(size(J1));
for i   = 1:length(V1)
    TXT(i)  = {sprintf('(%d, %d)',J1(i),J2(i))};
end

end

