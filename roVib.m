function [WAVELENGTH] = roVib(system, V1, V2, J1, J2, Fortrat)
% This function outputs the wavenumber of a certain ro-vibrational
% transitions between 2 different states of nitrogen given by (v1,J1) ->
% (v2, J2). 
%   ADDTIONAL NOTES ON USAGE:
%       INPUTS:-
%           - system:   System associated with transition between
%                           electronic states. Choose either of the
%                           following:
%                               1. '1+' First positive system.
%                               2. '2+' Second positive system.
%                               3. '1-' First negative system.
%           - V1:       Vibrational level of emitting state (integer array)
%           - V2:       Vibrational level of lower state (integer array)
%           - J1:       Rotational level of emitting state (integer array)
%           - J2:       Rotational level of lower state (integer array)
%           - FORTRAT:  Plot fortrat diagram if input is 1
%                           Issues: Will probably output an error if v1 and 
%                               v2 are arrays. Will need to look into it.
%       OUTPUTS:-   
%           - WAVELENGTH: Wavelength of the transitions of interest in nm.
%
%       If both V and J are arrays, they will need to be of the same
%       length. 


% Polynomimial coefficients derieved from R.R Laher's papers. The
% coefficients of each vibrational state have been organized as shown below.
%   Other notes: 
%       - Both Tv and Bv have units of cm^-1 and will have to be
%           converted to wavelength after the fact. 
%       - The polynomial coeffiencients for
%           rotational constants are only valid for v = 0-4.
%       - Te is the electronic energy of the vibrational state, Tv is the
%           vibrational energy of a specified vibrational level, and F is the 
%           rotational energy of a specified rotaitonal level. v here
%           corresponds to the net change in energy during transition. B is
%           the rotational constant
%       - Te1 and Be1 correspond to the upper state corresponding to
%           (v',J') states. Likewise, Te2 and Be2 correspond to the lower
%           states (v",J"). 
%
%         [A 3-Sig    , B 3-Pi    , C 3-Pi    , X 2-Sig   , B 2-Sig  ] 
T0      = [ 49754.8   ,  59306.8  ,  88977.9  ,  125667.5 ,  151233.5];
w       = [ 1460.48   ,  1734.38  ,  2047.18  ,  2207.37  ,  2420.83 ];
wx      = [ 13.775    ,  14.558   ,  28.445   ,  16.302   ,  23.851  ];
wy      = [-1.175e-2  ,  1.40e-2  ,  2.0883   , -2.67e-3  , -0.3587  ];
wz      = [ 1.41e-4   , -1.13e-3  , -5.350e-1 , -2.61e-3  , -6.192e-2];
wa      = [-7.29e-5   ,  0.0      ,  0.0      ,  3.7e-5   ,  0.0     ];

B0      = [ 1.45499   ,  1.63802  ,  1.8247   ,  1.93177  ,  2.0845   ];
a       = [ 1.8385e-2 ,  1.8302e-2,  1.868e-2 ,  1.900e-2 ,  2.132e-2 ];
ax      = [ 1.24e-5   , -8.4e-6   , -2.28e-3  , -1.91e-5  , -8.5e-4   ];
ay      = [-6.7e-6    , -3.4e-6   ,  7.33e-4  , -5.00e-6  ,  0        ];
az      = [ 0.0       ,  0.0      ,  0.0      ,  4.6e-8   ,  0.0      ];

switch system
    case '1+'
        % Constants for first positive system (corresponding to transtions 
        % from B 3-Gam to A 3-Sig electronic states)(i.e: 2->1)
        Te1     = T0(2) - w(2)/2 + wx(2)/4 - wy(2)/8 - wz(2)/16 - wa(2)/32;
        Te2     = T0(1) - w(1)/2 + wx(1)/4 - wy(1)/8 - wz(1)/16 - wa(1)/32;
        
        Tv1     = Te1 + w(2)  * (V1 + 0.5)    ...
                      - wx(2) * (V1 + 0.5).^2 ...
                      + wy(2) * (V1 + 0.5).^3 ...
                      + wz(2) * (V1 + 0.5).^4 ...
                      + wa(2) * (V1 + 0.5).^5;
        Tv2     = Te2 + w(1)  * (V2 + 0.5)    ...
                      - wx(1) * (V2 + 0.5).^2 ...
                      + wy(1) * (V2 + 0.5).^3 ...
                      + wz(1) * (V2 + 0.5).^4 ...
                      + wa(1) * (V2 + 0.5).^5;
                  
        Bv1     = B0(2) -  a(2) * (V1 + 0.5)    ...
                        + ax(2) * (V1 + 0.5).^2 ...
                        + ay(2) * (V1 + 0.5).^3 ... 
                        + az(2) * (V1 + 0.5).^4; 
        Bv2     = B0(1) -  a(1) * (V2 + 0.5)    ...
                        + ax(1) * (V2 + 0.5).^2 ...
                        + ay(1) * (V2 + 0.5).^3 ...
                        + az(1) * (V2 + 0.5).^4;

    case '2+'
        % Constants for second positive system (corresponding to transtions 
        % from C 3-Gam to B 3-Gam electronic states). (i.e: 3->2)
        Te1     = T0(3) - w(3)/2 + wx(3)/4 - wy(3)/8 - wz(3)/16 - wa(3)/32;
        Te2     = T0(2) - w(2)/2 + wx(2)/4 - wy(2)/8 - wz(2)/16 - wa(2)/32;
        
        Tv1     = Te1 + w(3)  * (V1 + 0.5)    ...
                      - wx(3) * (V1 + 0.5).^2 ...
                      + wy(3) * (V1 + 0.5).^3 ...
                      + wz(3) * (V1 + 0.5).^4 ...
                      + wa(3) * (V1 + 0.5).^5;
        Tv2     = Te2 + w(2)  * (V2 + 0.5)    ...
                      - wx(2) * (V2 + 0.5).^2 ...
                      + wy(2) * (V2 + 0.5).^3 ...
                      + wz(2) * (V2 + 0.5).^4 ...
                      + wa(2) * (V2 + 0.5).^5;
                  
        Bv1     = B0(3) -  a(3) * (V1 + 0.5)    ...
                        + ax(3) * (V1 + 0.5).^2 ...
                        + ay(3) * (V1 + 0.5).^3 ...
                        + az(3) * (V1 + 0.5).^4; 
        Bv2     = B0(2) -  a(2) * (V2 + 0.5)    ...
                        + ax(2) * (V2 + 0.5).^2 ...
                        + ay(2) * (V2 + 0.5).^3 ...
                        + az(2) * (V2 + 0.5).^4;

    case '1-'
        % Constants for second positive system (corresponding to transtions 
        % from B 2-Sig to X 2-Sig electronic states). (i.e: 5->4)
        Te1     = T0(5) - w(5)/2 + wx(5)/4 - wy(5)/8 - wz(5)/16 - wa(5)/32;
        Te2     = T0(4) - w(4)/2 + wx(4)/4 - wy(4)/8 - wz(4)/16 - wa(4)/32;
        
        Tv1     = Te1 + w(5)  * (V1 + 0.5)    ...
                      - wx(5) * (V1 + 0.5).^2 ...
                      + wy(5) * (V1 + 0.5).^3 ...
                      + wz(5) * (V1 + 0.5).^4 ...
                      + wa(5) * (V1 + 0.5).^5;
        Tv2     = Te2 + w(4)  * (V2 + 0.5)    ...
                      - wx(4) * (V2 + 0.5).^2 ...
                      + wy(4) * (V2 + 0.5).^3 ...
                      + wz(4) * (V2 + 0.5).^4 ...
                      + wa(4) * (V2 + 0.5).^5;
                  
        Bv1     = B0(5) -  a(5) * (V1 + 0.5)    ...
                        + ax(5) * (V1 + 0.5).^2 ...
                        + ay(5) * (V1 + 0.5).^3 ...
                        + az(5) * (V1 + 0.5).^4; 
        Bv2     = B0(4) -  a(4) * (V2 + 0.5)    ...
                        + ax(4) * (V2 + 0.5).^2 ...
                        + ay(4) * (V2 + 0.5).^3 ...
                        + az(4) * (V2 + 0.5).^4;

    otherwise 
        error('Incorrect transistion requested or input is not a string');
end

F1      = Bv1 .* J1.*(J1 + 1);
F2      = Bv2 .* J2.*(J2 + 1);
v       = (Tv1 - Tv2) + (F1 - F2);


WAVELENGTH   = 1e7./v;

% Plot the fortrat diagram for the first 50 rotational levels for user
%   defined transition from (v', v"). Divide up the data into two sets. One
%   for the P-branch transitions (i.e. transition with delJ = +1) and
%   R-branch transititions (i.e. transitions with delJ = -1). 
%       Other notes: 
%           -j2 is this lower rotational state while j1 is the upper
%           rotational state.  k is the wavenumber (in cm^-1). wl is the
%           wavelength (in nm). 
if Fortrat  == 1
    wl      = zeros(50,1);
    %P-branch transitions
    for i = 1:50
        j2      = i;
        j1      = j2-1;
        f1      = Bv1 * j1 * (j1 + 1);
        f2      = Bv2 * j2 * (j2 + 1);
        k       = (Tv1 - Tv2) + (f1 - f2);
        wl(i)   = 1e7/k;
    end 
    scatter (wl, linspace(1,50,50),'r','filled','DisplayName','P-Branch'); 
    hold on;
    
    %R-branch transitions
    for i = 1:50
        j2      = i-1;
        j1      = j2+1;
        f1      = Bv1 * j1 * (j1 + 1);
        f2      = Bv2 * j2 * (j2 + 1);
        k       = (Tv1 - Tv2) + (f1 - f2);
        wl(i)   = 1e7/k;
    end
    
    scatter (wl, linspace(1,50,50),'k','filled','DisplayName','R-Branch');
    
    xlabel('Wavelength (nm)');
    ylabel('Rotational quantum number, J''''');
    title(sprintf('Fortrat Diagram for (v'', v'''') = (%d,%d)',V1, V2));
    legend('show')
end

    