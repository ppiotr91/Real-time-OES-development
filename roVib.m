function [WAVELENGTH, BV1, BV2] = roVib(system, V1, V2, J1, J2)
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
%
%       OUTPUTS:-   
%           - WAVELENGTH: Wavelength of the transitions of interest in nm.
%           - BV1:        Rotational constant 
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
%           rotational energy of a specified rotational level. v here
%           corresponds to the net change in energy during transition. B is
%           the rotational constant
%       - Te1 and Be1 correspond to the upper state corresponding to
%           (v',J') states. Likewise, Te2 and Be2 correspond to the lower
%           states (v",J"). 
%
%         [ A 3-Sig    , B 3-Pi    , C 3-Pi    , X 2-Sig   , B 2-Sig  ] 
T0      = [  50203.63  ,  59619.35 ,  89136.88 ,  0.0      ,  25461.46];
w       = [  1460.64   ,  1733.39  ,  2047.178 ,  2207.00  ,  2419.84 ];
wx      = [  13.872    ,  14.122   ,  28.445   ,  16.10    ,  23.189  ];
wy      = [  0.0103    , -0.0569   ,  2.08833  , -0.040    , -0.5375  ];
wz      = [ -0.00197   ,  0.00361  , -0.5350   ,  0.0      , -0.0495  ];
wa      = [  0.0       ,  0.0      ,  0.0      ,  0.0      ,  0.0     ];

B0      = [  1.4546    ,  1.63745  ,  1.82473  ,  1.93176  ,  2.0845   ];
a       = [  0.0180    ,  0.01791  ,  1.868e-2 ,  0.01881  ,  0.024    ];
ax      = [ -3.3e-5    , -7.7e-5   , -2.28e-3  ,  0.0      ,  0.0      ];
ay      = [  0.0       ,  0.0      ,  7.33e-4  ,  0.0      ,  0.0      ];
az      = [  0.0       ,  0.0      , -1.5e-04  ,  0.0      ,  0.0      ];

De      = [  5.76e-6   , 5.9e-6    ,  0.0      ,  6.10e-6  ,  6.17e-6  ];

switch system
    case '1+'
        % Constants for first positive system (corresponding to transtions 
        % from B 3-Gam to A 3-Sig electronic states)(i.e: 2->1)
        Te1     = T0(2);
        Te2     = T0(1);
        
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
                    
        Dv1     = De(2);
        Dv2     = De(1);

    case '2+'
        % Constants for second positive system (corresponding to transtions 
        % from C 3-Gam to B 3-Gam electronic states). (i.e: 3->2)
        Te1     = T0(3);
        Te2     = T0(2) - 8.3877 - 3.1487 + 3.2543;
        
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
                        + az(2) * (V2 + 0.5).^4 ...
                        - 6.95e-4;
        Dv1     = De(3);
        Dv2     = De(2);

    case '1-'
        % Constants for second positive system (corresponding to transtions 
        % from B 2-Sig to X 2-Sig electronic states). (i.e: 5->4)
        Te1     = T0(5);
        Te2     = T0(4);
        
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
                    
        Dv1     = De(5);
        Dv2     = De(4);

    otherwise 
        error('Incorrect transistion requested or input is not a string');
end
F1      = Bv1 .* J1.*(J1 + 1) - Dv1 * J1.^2 .* (J1 + 1).^2;
F2      = Bv2 .* J2.*(J2 + 1) - Dv2 * J2.^2 .* (J2 + 1).^2;
v       = (Tv1 - Tv2) + (F1 - F2);

switch nargout
    case 1
        WAVELENGTH   = 1e7./v;
    case 2
        WAVELENGTH   = 1e7./v;
        BV1 = Bv1;
    case 3
        WAVELENGTH   = 1e7./v;
        BV1 = Bv1;
        BV2 = Bv2;
end

end

    