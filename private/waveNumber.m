function [data] = waveNumber(syst, V1, V2, Jf, constants)
con1 	= constants.con1(V1+1,:);
con2 	= constants.con2(V2+1,:);

Te1     = con1(1);
Gv1     = con1(2);
Bv1     = con1(3);
Dv1     = con1(4);
Y1      = con1(5);
g1      = con1(6);
epsi1 	= con1(7);

Te2     = con2(1);
Gv2     = con2(2);
Bv2     = con2(3);
Dv2     = con2(4);
Y2      = con2(5);
g2      = con2(6);
epsi2 	= con2(7); 

EnerDiff    = (Te1 - Te2) + (Gv1 - Gv2);

switch syst
    case {'2+', 'swan'}
        J1(:,1)     = (0:Jf)';
        J2(:,1)     = (0:Jf)';
        
        y1(1,1)     = Y1*(Y1 - 4) + 4/3;
        y1(2,1)     = Y1*(Y1 - 1) - 4/9;
        y2(1,1)     = Y2*(Y2 - 4) + 4/3;
        y2(2,1)     = Y2*(Y2 - 1) - 4/9;
        jj1         = J1.*(J1 + 1);
        jj2         = J2.*(J2 + 1);
        z1(:,1)     = y1(1) + 4*jj1;
        z1(:,2)     = 2/3 * (y1(2) - 2*jj1)./z1(:,1);
        z2(:,1)     = y2(1) + 4*jj2;
        z2(:,2)     = 2/3 * (y2(2) - 2*jj2)./z2(:,1);
        F01(:,1)    = Bv1 * jj1;
        F02(:,1)    = Bv2 * jj2;
        
        %F  = [f1, f2, f3] where f is a substate
        F1(:,1)     = F01 + Bv1 * (-z1(:,1).^0.5 - z1(:,2)) - Dv1*(J1 - 0.5).^4;
        F1(:,2)     = F01 + Bv1 * (2 * z1(:,2))             - Dv1*(J1 + 0.5).^4;
        F1(:,3)     = F01 + Bv1 * ( z1(:,1).^0.5 - z1(:,2)) - Dv1*(J1 + 1.5).^4;
        
        F2(:,1)     = F02 + Bv2 * (-z2(:,1).^0.5 - z2(:,2)) - Dv2*(J2 - 0.5).^4;
        F2(:,2)     = F02 + Bv2 * (2 * z2(:,2))             - Dv2*(J2 + 0.5).^4;
        F2(:,3)     = F02 + Bv2 * ( z2(:,1).^0.5 - z2(:,2)) - Dv2*(J2 + 1.5).^4;
        
        range = defineRange(syst,J1);
        
	case 'nhA-X'
		J1(:,1)     = (0:Jf)';
		y1(1,1)     = Y1*(Y1 - 4) + 4/3;
        y1(2,1)     = Y1*(Y1 - 1) - 4/9;
		jj1         = J1.*(J1 + 1);
		z1(:,1)     = y1(1) + 4*jj1;
		z1(:,2)     = 2/3 * (y1(2) - 2*jj1)./z1(:,1);
		F01(:,1)    = Bv1 * jj1;
		
		F1(:,1)     = F01 + Bv1 * (-z1(:,1).^0.5 - z1(:,2)) - Dv1*(J1 - 0.5).^4;
        F1(:,2)     = F01 + Bv1 * (2 * z1(:,2))             - Dv1*(J1 + 0.5).^4;
        F1(:,3)     = F01 + Bv1 * ( z1(:,1).^0.5 - z1(:,2)) - Dv1*(J1 + 1.5).^4;
		
		N2  = (0:Jf)';
		J2 	= [N2-1, N2 ,N2+1];
		nn2 = N2.*(N2+1);
		f0  = Bv2*nn2 - Dv2*nn2.^2;
		
		w2(:,1) 	=    Bv2*(2*N2+3) - epsi2 - ( Bv2^2*(2*N2+3).^2 + epsi2^2 - 2*epsi2*Bv2 ).^0.5 + g2*(N2+1);
		w2(:,2) 	= 0;
		w2(:,3) 	= -1*Bv2*(2*N2-1) - epsi2 + ( Bv2^2*(2*N2-1).^2 + epsi2^2 - 2*epsi2*Bv2 ).^0.5 - g2*N2;
		F2 			= f0 + w2;
		
		range 	= defineRange(syst,J1);
		
		data.N2 = N2;
		
    case {'1-', 'CNvio'}
        N1  = (0:Jf)';
        N2  = (0:Jf)';
        J1 	= [N1+0.5, N1-0.5];
        J2  = [N2+0.5, N2-0.5];
        nn1 = N1.*(N1+1);
        nn2 = N2.*(N2+1);
        
        %F  = [f1, f2] where f is a substate
        F1(:,1) = Bv1*nn1 - Dv1*nn1.^2 + 0.5*g1.*N1;
        F2(:,1) = Bv2*nn2 - Dv1*nn2.^2 + 0.5*g2.*N2;
        F1(:,2) = Bv1*nn1 - Dv1*nn1.^2 - 0.5*g1.*(N1+1);
        F2(:,2) = Bv2*nn2 - Dv2*nn2.^2 - 0.5*g2.*(N2+1);
        
        range = defineRange(syst,N1);

        data.N1     = N1;
        data.N2     = N2; 
		
    otherwise
        error('Incorrect transistion requested or input is not a string');
end

p = calculateWaveNumbers(F1 + EnerDiff, F2, range.p.j1, range.p.j2);
q = calculateWaveNumbers(F1 + EnerDiff, F2, range.q.j1, range.q.j2);
r = calculateWaveNumbers(F1 + EnerDiff, F2, range.r.j1, range.r.j2); 
data.wlP    = wNum2wLength(p);
data.wlQ    = wNum2wLength(q);
data.wlR    = wNum2wLength(r);

data.syst   = syst;

data.Te1    = Te1;
data.Gv1    = Gv1;
data.Bv1    = Bv1;
data.Dv1    = Dv1;
data.Y1     = Y1;
data.g1     = g1;

data.Te2    = Te2;
data.Gv2    = Gv2;
data.Bv2    = Bv2;
data.Dv2    = Dv2;
data.Y2     = Y2;
data.g2     = g2;

data.V1     = V1;
data.V2     = V2;
data.J1     = J1;
data.J2     = J2;

data.F1     = F1;
data.F2     = F2;

data.range  = range;

data.wlbh   = findBandHead(data.wlP);

 
    function [wavelength] = wNum2wLength(v)
        lambda      = 1e4./v;             %wavelength [um]
        sigma       = 1./lambda;          % [um^-1]
        sigma       = sigma.^2;
        nAir        = (8342.13 + 2406030./(130-sigma) + 15997./(38.9 - sigma)) * 1e-8 + 1;
        wavelength(:,:,:)       = 1e7./v.* 1./nAir;
        wavelength(isinf(wavelength) == 1) = 0; 
    end

    function [range] = defineRange(syst, J1)
        %define the indicies of the rotational quntum number, J, for both
        %upper and lower states. Only define values for transitions of
        %inter
        switch syst
            case '2+'
                lenf = length(J1);
                range.p.j1  = zeros(lenf,3,3);
                range.p.j2  = zeros(lenf,3,3);
                range.q.j1  = zeros(lenf,3,3);
                range.q.j2  = zeros(lenf,3,3);
                range.r.j1  = zeros(lenf,3,3);
                range.r.j2  = zeros(lenf,3,3);
                
                %P-Branch:
                lenf = length(J1);
                range.p.j1(1:lenf-1,1,1)    = 1:lenf-1;
                range.p.j1(2:lenf-1,2,2)    = 2:lenf-1;
                range.p.j1(3:lenf-1,3,3)    = 3:lenf-1;
                range.p.j2(2:lenf,1,1)      = 2:lenf;
                range.p.j2(3:lenf,2,2)      = 3:lenf;
                range.p.j2(4:lenf,3,3)      = 4:lenf;
                %Q-Branch
                range.q.j1(2:lenf,2,2)  = 2:lenf;
                range.q.j1(3:lenf,3,3)  = 3:lenf;
                range.q.j2              = range.q.j1;
                %R-Branch
                range.r.j1(2:lenf,1,1)      = 2:lenf;
                range.r.j1(3:lenf,2,2)      = 3:lenf;
                range.r.j1(4:lenf,3,3)      = 4:lenf;
                range.r.j2(1:lenf-1,1,1)    = 1:lenf-1;
                range.r.j2(2:lenf-1,2,2)    = 2:lenf-1;
                range.r.j2(3:lenf-1,3,3)    = 3:lenf-1;
				
            case 'swan' 
                lenf = length(J1);
                range.p.j1  = zeros(lenf,3,3);
                range.p.j2  = zeros(lenf,3,3);
                range.q.j1  = zeros(lenf,3,3);
                range.q.j2  = zeros(lenf,3,3);
                range.r.j1  = zeros(lenf,3,3);
                range.r.j2  = zeros(lenf,3,3);
                
                %P-Branch
                range.p.j1(3:lenf-1,1,1)    = 3:lenf-1;    
                range.p.j1(2:lenf-1,2,2)    = 2:lenf-1;    
                range.p.j1(1:lenf-1,3,3)    = 1:lenf-1;    
                range.p.j2(4:lenf,1,1)      = 4:lenf;
                range.p.j2(3:lenf,2,2)      = 3:lenf;
                range.p.j2(2:lenf,3,3)      = 2:lenf;
                %Q-Branch
                range.q.j1(3:lenf,1,1)      = 3:lenf;
                range.q.j1(2:lenf,2,2)      = 2:lenf;      
                range.q.j2                  = range.q.j1;
                %R-Branch
                range.r.j1(4:lenf,1,1)      = 4:lenf;
                range.r.j1(3:lenf,2,2)      = 3:lenf;
                range.r.j1(2:lenf,3,3)      = 2:lenf;
                range.r.j2(3:lenf-1,1,1)    = 3:lenf-1;
                range.r.j2(2:lenf-1,2,2)    = 2:lenf-1;
                range.r.j2(1:lenf-1,3,3)    = 1:lenf-1;
				
			case 'chB-X'
				lenf = length(J1);
                range.p.j1  = zeros(lenf,3,3);
                range.p.j2  = zeros(lenf,3,3);
                range.q.j1  = zeros(lenf,3,3);
                range.q.j2  = zeros(lenf,3,3);
                range.r.j1  = zeros(lenf,3,3);
                range.r.j2  = zeros(lenf,3,3);
                
                %P-Branch
                range.p.j1(3:lenf-1,1,1)    = 3:lenf-1;    
                range.p.j1(2:lenf-1,2,2)    = 2:lenf-1;    
                range.p.j1(1:lenf-1,3,3)    = 1:lenf-1;    
                range.p.j2(4:lenf,1,1)      = 4:lenf;
                range.p.j2(3:lenf,2,2)      = 3:lenf;
                range.p.j2(2:lenf,3,3)      = 2:lenf;
                %Q-Branch
                range.q.j1(3:lenf,1,1)      = 3:lenf;
                range.q.j1(2:lenf,2,2)      = 2:lenf;      
                range.q.j2                  = range.q.j1;
                %R-Branch
                range.r.j1(4:lenf,1,1)      = 4:lenf;
                range.r.j1(3:lenf,2,2)      = 3:lenf;
                range.r.j1(2:lenf,3,3)      = 2:lenf;
                range.r.j2(3:lenf-1,1,1)    = 3:lenf-1;
                range.r.j2(2:lenf-1,2,2)    = 2:lenf-1;
                range.r.j2(1:lenf-1,3,3)    = 1:lenf-1;
				
            case {'1-', 'CNvio'}
                %j1 and j2 will be used as stand-ins for n1 and n2 for the
                %ranged function!
                lenf = length(J1); 
                range.p.j1  = zeros(lenf,2,2);
                range.p.j2  = zeros(lenf,2,2);
                range.q.j1  = zeros(lenf,2,2);
                range.q.j2  = zeros(lenf,2,2);
                range.r.j1  = zeros(lenf,2,2);
                range.r.j2  = zeros(lenf,2,2);
                
                %P-Branch
                range.p.j1(1:lenf-1,1,1)    = 1:lenf-1;
                range.p.j1(2:lenf-1,2,2)    = 2:lenf-1;
                range.p.j2(2:lenf,1,1)      = 2:lenf;
                range.p.j2(3:lenf,2,2)      = 3:lenf;
                %Q-Branch
                range.q.j1(2:lenf,2,1)  = 2:lenf;
                range.q.j1(2:lenf,1,2)  = 2:lenf;
                range.q.j2              = range.q.j1;
                %R-Branch
                range.r.j1(2:lenf,1,1)      = 2:lenf;
                range.r.j1(3:lenf,2,2)      = 3:lenf;
                range.r.j2(1:lenf-1,1,1)    = 1:lenf-1;
                range.r.j2(2:lenf-1,2,2)    = 2:lenf-1;
        end
    end
    
    function [out] = calculateWaveNumbers(F1,F2,range1,range2)
        %I denotes upper level
        %J denotes lower level
        maxI    = size(F1,2);
        maxJ    = size(F2,2);
        out     = zeros(size(F1(:,1),1), maxJ, maxI);
        
        for i   = 1:maxI
            for j = 1:maxJ
                %sample equation with range1 and range2 values:
                %   x1(1:end-1,1)     = F1(1:end-1,1) - F2(2:end,1) + EnerDiff;
                if  sum(range1(:,j,i)) == 0; continue; end
                if  sum(range2(:,j,i)) == 0; continue; end
                r1  = range1(range1(:,j,i) ~= 0 , j,i);
                r2  = range2(range2(:,j,i) ~= 0 , j,i);
                out(r1,j,i)     = F1(r1,i) - F2(r2,j);
            end   
        end
    end  
    
    function [out] = findBandHead(wl)
        out     = max(max(max(wl)));
    end
end
