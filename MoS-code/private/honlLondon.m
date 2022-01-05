function [Sp,Sq,Sr] = honlLondon(bandData)
system  = bandData.syst;
J1      = bandData.J1; J1   = [J1; J1(end,:) + 1];
J2      = bandData.J2; J2   = [J2; J2(end,:) + 1];
Y1      = bandData.Y1;
Y2      = bandData.Y2;

switch system
    case {'2+','swan'}
        lambda  = 1;
        Sp      = zeros(length(J1),3,3);
        Sq      = zeros(length(J1),3,3);
        Sr      = zeros(length(J1),3,3);  
        
        const 	= hlTripletConst('inter', lambda, Y1, Y2, J1, J2);
        
        %P-Branch
        range1  = bandData.range.p.j1;
        range2  = bandData.range.p.j2;
        for i   = 1:3
            for j   = 1:3
                if  sum(range1(:,j,i)) == 0; continue; end
                if  sum(range2(:,j,i)) == 0; continue; end
                r1  = range1(range1(:,j,i) ~= 0 , j,i); rr1 = [r1; max(r1) + 1];
                r2  = range2(range2(:,j,i) ~= 0 , j,i); rr2 = [r2; max(r2) + 1];
                branch  = sprintf('p%d%d',i,j);
                Sp(rr1,j,i)  = hlTripletCalc(branch,const, lambda, Y1, Y2, J1, J2, rr1, rr2);
            end
        end
        %Q-Branch
        range1  = bandData.range.q.j1;
        range2  = bandData.range.q.j2;
        for i   = 1:3
            for j   = 1:3
                if  sum(range1(:,j,i)) == 0; continue; end
                if  sum(range2(:,j,i)) == 0; continue; end
                r1  = range1(range1(:,j,i) ~= 0 ,j,i); rr1 = [r1; max(r1) + 1];
                r2  = range2(range2(:,j,i) ~= 0 ,j,i); rr2 = [r2; max(r2) + 1];
                branch  = sprintf('q%d%d',i,j);
                Sq(rr1,j,i)  = hlTripletCalc(branch,const, lambda, Y1, Y2, J1, J2, rr1, rr2);
            end
        end
        %R-Branch
        range1  = bandData.range.r.j1;
        range2  = bandData.range.r.j2;
        for i   = 1:3
            for j   = 1:3
                if  sum(range1(:,j,i)) == 0; continue; end
                if  sum(range2(:,j,i)) == 0; continue; end
                r1  = range1(range1(:,j,i) ~= 0 ,j,i); rr1 = [r1; max(r1) + 1];
                r2  = range2(range2(:,j,i) ~= 0 ,j,i); rr2 = [r2; max(r2) + 1];
                branch  = sprintf('r%d%d',i,j);
                Sr(rr1,j,i)  = hlTripletCalc(branch,const, lambda, Y1, Y2, J1, J2, rr1, rr2);
            end
        end
	
	case 'nhA-X'
		lambda 	= 0; %Smaller of the 2 lambda is used for transitions between different types of electronic states
        Sp      = zeros(length(J1),3,3);
        Sq      = zeros(length(J1),3,3);
        Sr      = zeros(length(J1),3,3);  
		
		const 	= hlTripletConst('inter', lambda, Y1, Y2, J1, J2);

		%P-Branch
        range1  = bandData.range.p.j1;
        range2  = bandData.range.p.j2;
        for i   = 1:3
            for j   = 1:3
                if  sum(range1(:,j,i)) == 0; continue; end
                if  sum(range2(:,j,i)) == 0; continue; end
                r1  = range1(range1(:,j,i) ~= 0 , j,i); rr1 = [r1; max(r1) + 1];
                r2  = range2(range2(:,j,i) ~= 0 , j,i); rr2 = [r2; max(r2) + 1];
                branch  = sprintf('p%d%d',i,j);
				
                j2          = J2(:,j);   
                Sp(rr1,j,i)  = hlTripletCalcDel1(branch,const, lambda, Y1, Y2, J1, j2, rr1, rr2);
            end
        end
		%Q-Branch
        range1  = bandData.range.q.j1;
        range2  = bandData.range.q.j2;
        for i   = 1:3
            for j   = 1:3
                if  sum(range1(:,j,i)) == 0; continue; end
                if  sum(range2(:,j,i)) == 0; continue; end
                r1  = range1(range1(:,j,i) ~= 0 , j,i); rr1 = [r1; max(r1) + 1];
                r2  = range2(range2(:,j,i) ~= 0 , j,i); rr2 = [r2; max(r2) + 1];
                branch  = sprintf('q%d%d',i,j);
				
                j2          = J2(:,j);   
                Sq(rr1,j,i)  = hlTripletCalcDel1(branch,const, lambda, Y1, Y2, J1, j2, rr1, rr2);
            end
        end
		%R-Branch
        range1  = bandData.range.r.j1;
        range2  = bandData.range.r.j2;
        for i   = 1:3
            for j   = 1:3
                if  sum(range1(:,j,i)) == 0; continue; end
                if  sum(range2(:,j,i)) == 0; continue; end
                r1  = range1(range1(:,j,i) ~= 0 , j,i); rr1 = [r1; max(r1) + 1];
                r2  = range2(range2(:,j,i) ~= 0 , j,i); rr2 = [r2; max(r2) + 1];
                branch  = sprintf('r%d%d',i,j);
				
                j2          = J2(:,j);   
                Sr(rr1,j,i)  = hlTripletCalcDel1(branch,const, lambda, Y1, Y2, J1, j2, rr1, rr2);
            end
        end
	
    case '1-'
        lambda  = 0;
        Sp      = zeros(length(J1),2,2);
        Sq      = zeros(length(J1),2,2);
        Sr      = zeros(length(J1),2,2);
        
        %P-Branch
        range1  = bandData.range.p.j1;
        range2  = bandData.range.p.j2;
        for i   = 1:2
            for j   = 1:2
                if  sum(range1(:,j,i)) == 0; continue; end
                if  sum(range2(:,j,i)) == 0; continue; end
                r1  = range1(range1(:,j,i) ~= 0 ,j,i); rr1 = [r1; max(r1) + 1];
                r2  = range2(range2(:,j,i) ~= 0 ,j,i); rr2 = [r2; max(r2) + 1];
                
                branch      = sprintf('p%d%d',i,j);
                j1          = J1(:,i);
                j2          = J2(:,j);                
                const       = hlDoubletConst('b', lambda, Y1,  Y2, j1, j2);
                Sp(rr1,j,i) = hlDoubletCalc(branch,const, lambda, j1, j2, rr1, rr2);
            end
        end
        %Q-Branch
        range1  = bandData.range.q.j1;
        range2  = bandData.range.q.j2;
        for i   = 1:2
            for j   = 1:2
                if  sum(range1(:,j,i)) == 0; continue; end
                if  sum(range2(:,j,i)) == 0; continue; end
                r1  = range1(range1(:,j,i) ~= 0 ,j,i); rr1 = [r1; max(r1) + 1];
                r2  = range2(range2(:,j,i) ~= 0 ,j,i); rr2 = [r2; max(r2) + 1];
                
                branch      = sprintf('q%d%d',i,j);
                j1          = J1(:,i);
                j2          = J2(:,j);                
                const       = hlDoubletConst('b', lambda, Y1,  Y2, j1, j2);
                Sq(rr1,j,i) = hlDoubletCalc(branch,const, lambda, j1, j2, rr1, rr2);
            end
        end
        %R-Branch
        range1  = bandData.range.r.j1;
        range2  = bandData.range.r.j2;
        for i   = 1:2
            for j   = 1:2
                if  sum(range1(:,j,i)) == 0; continue; end
                if  sum(range2(:,j,i)) == 0; continue; end
                r1  = range1(range1(:,j,i) ~= 0 ,j,i); rr1 = [r1; max(r1) + 1];
                r2  = range2(range2(:,j,i) ~= 0 ,j,i); rr2 = [r2; max(r2) + 1];
                
                branch      = sprintf('r%d%d',i,j);
                j1          = J1(:,i);
                j2          = J2(:,j);                
                const       = hlDoubletConst('b', lambda, Y1,  Y2, j1, j2);
                Sr(rr1,j,i) = hlDoubletCalc(branch,const, lambda, j1, j2, rr1, rr2);
            end
        end
                
        %alternation correction (see Hanson)
        alter   = ones(size(Sp));
        alter(mod(bandData.N1,2) == 1,:,:) = 2/3; %odd N'
        alter(mod(bandData.N1,2) == 0,:,:) = 1/3; %even N"
        
        Sp  = Sp.*alter;
        Sq  = Sq.*alter;
        Sr  = Sr.*alter;
                
    case 'CNvio'
        lambda  = 0;
        Sp      = zeros(length(J1),2,2);
        Sq      = zeros(length(J1),2,2);
        Sr      = zeros(length(J1),2,2);
        
        %P-Branch
        range1  = bandData.range.p.j1;
        range2  = bandData.range.p.j2;
        for i   = 1:2
            for j   = 1:2
                if  sum(range1(:,j,i)) == 0; continue; end
                if  sum(range2(:,j,i)) == 0; continue; end
                r1  = range1(range1(:,j,i) ~= 0 ,j,i); rr1 = [r1; max(r1) + 1];
                r2  = range2(range2(:,j,i) ~= 0 ,j,i); rr2 = [r2; max(r2) + 1];
                
                branch      = sprintf('p%d%d',i,j);
                j1          = J1(:,i);
                j2          = J2(:,j);                
                const       = hlDoubletConst('b', lambda, Y1,  Y2, j1, j2);
                Sp(rr1,j,i) = hlDoubletCalc(branch,const, lambda, j1, j2, rr1, rr2);
            end
        end
        %Q-Branch
        range1  = bandData.range.q.j1;
        range2  = bandData.range.q.j2;
        for i   = 1:2
            for j   = 1:2
                if  sum(range1(:,j,i)) == 0; continue; end
                if  sum(range2(:,j,i)) == 0; continue; end
                r1  = range1(range1(:,j,i) ~= 0 ,j,i); rr1 = [r1; max(r1) + 1];
                r2  = range2(range2(:,j,i) ~= 0 ,j,i); rr2 = [r2; max(r2) + 1];
                
                branch      = sprintf('q%d%d',i,j);
                j1          = J1(:,i);
                j2          = J2(:,j);                
                const       = hlDoubletConst('b', lambda, Y1,  Y2, j1, j2);
                Sq(rr1,j,i) = hlDoubletCalc(branch,const, lambda, j1, j2, rr1, rr2);
            end
        end
        %R-Branch
        range1  = bandData.range.r.j1;
        range2  = bandData.range.r.j2;
        for i   = 1:2
            for j   = 1:2
                if  sum(range1(:,j,i)) == 0; continue; end
                if  sum(range2(:,j,i)) == 0; continue; end
                r1  = range1(range1(:,j,i) ~= 0 ,j,i); rr1 = [r1; max(r1) + 1];
                r2  = range2(range2(:,j,i) ~= 0 ,j,i); rr2 = [r2; max(r2) + 1];
                
                branch      = sprintf('r%d%d',i,j);
                j1          = J1(:,i);
                j2          = J2(:,j);                
                const       = hlDoubletConst('b', lambda, Y1,  Y2, j1, j2);
                Sr(rr1,j,i) = hlDoubletCalc(branch,const, lambda, j1, j2, rr1, rr2);
            end
        end
end

Sp(isnan(Sp) == 1 | isinf(Sp) == 1)  = 0;
Sq(isnan(Sq) == 1 | isinf(Sq) == 1)  = 0;
Sr(isnan(Sr) == 1 | isinf(Sr) == 1)  = 0;

Sp(end,:,:) = [];
Sq(end,:,:) = [];
Sr(end,:,:) = [];


    function [HLconst] 	= hlDoubletConst(type, lambda, Y1, Y2, J1, J2)
        J 	= [J1(:,1), J2(:,1)];
        I 	= ones(size(J1));
        Y 	= [Y1.*I, Y2.*I];
        
        switch type
            case 'a'
                HLconst.up 	= 2*lambda*(Y-2);
                HLconst.un 	= zeros(size(J));
                HLconst.cp 	= 2*lambda^2*(Y-2).^2;
                HLconst.cn 	= 2*(J-lambda+0.5).*(J+lambda+0.5);
                
            case 'b'
                HLconst.up 	= 2*(J-lambda+0.5);
                HLconst.un 	= 2*(J+lambda+0.5);
                HLconst.cp 	= 4*(J+0.5).*(J-lambda+0.5);
                HLconst.cn 	= 4*(J+0.5).*(J+lambda+0.5);
        end
        
    end

    function [HLconst]	= hlTripletConst(type, lambda, Y1, Y2, J1, J2)
        J 	= [J1, J2];
        I 	= ones(size(J1));
        
        switch type
            case 'inter'
                m(:,1) 	= lambda^2*Y1*(Y1-4)*I;
                m(:,2)  = lambda^2*Y2*(Y2-4)*I;
                n(:,1)  = lambda*(Y1-2)*I;
                n(:,2)  = lambda*(Y2-2)*I;
                
                HLconst.u1p = (m + 4*J.^2).^0.5 +  n;
                HLconst.u1n = (m + 4*J.^2).^0.5 -  n;
                HLconst.u3p = (m + 4*(J+1).^2).^0.5 + n;
                HLconst.u3n = (m + 4*(J+1).^2).^0.5 - n;
                
				if (Y1 && Y2 > 0)
					HLconst.c1 = m.*(J - lambda + 1).*(J + lambda) + 2*(2*J + 1).*(J - lambda).*J.*(J + lambda);
					HLconst.c2 = m + 4.*J.*(J+lambda);
					HLconst.c3 = m.*(J + lambda + 1).*(J - lambda) + 2*(2*J + 1).*(J - lambda + 1).*(J + 1).*(J + lambda + 1);
                elseif (Y1 && Y2 <0)
					HLconst.c2 = m + 4.*J.*(J+lambda);
					lambda  = -lambda;
					HLconst.c1 = m.*(J - lambda + 1).*(J + lambda) + 2*(2*J + 1).*(J - lambda).*J.*(J + lambda);
					HLconst.c3 = m.*(J + lambda + 1).*(J - lambda) + 2*(2*J + 1).*(J - lambda + 1).*(J + 1).*(J + lambda + 1);
				else 
					error('honlLondon\hlTriplet: Y1 and Y2 have opposite signs! Require further evalution for this case...')
				end
				
            case 'a'
                m(:,1) 	= lambda^2*(Y1-2)^2 * I;
                m(:,2) 	= lambda^2*(Y2-2)^2 * I;
                n(:,1) 	= lambda*(Y1-2) * I;
                n(:,2) 	= lambda*(Y2-2) * I;
                
                HLconst.u1p = 2*n;
                HLconst.u1n = zeros(size(J));
                HLconst.u3p = HLconst.u1p;
                HLconst.u3n = HLconst.u1n;
                
                HLconst.c1 	= m.*(J-lambda+1).*(J+lambda);
                HLconst.c2 	= m;
                HLconst.c3 	= m.*(J+lambda+1).*(J-lambda);
                
            case 'b'
                HLconst.u1p = 2*(J-lambda);
                HLconst.u1n = 2*(J+lambda);
                HLconst.u3p = 2*(J-lambda+1);
                HLconst.u3n = 2*(J+lambda+1);
                
                HLconst.c1 	= 2*(2*J+1).*(J-lambda).*J.*(J+lambda);
                HLconst.c2 	= 4*J.*(J+1);
                HLconst.c3 	= 2*(2*J+1).*(J-lambda+1).*(J+1).*(J+lambda+1);
        end
    end
	
    function [out] = hlDoubletCalc(branch, const, lambda, J1, J2, range1, range2)
        J1 		= J1(range1,1);
        J2 		= J2(range2,1);
        uPos    = [const.up(range1,1), const.up(range2,2)];
        uNeg    = [const.un(range1,1), const.un(range2,2)];
        cPos    = [const.cp(range1,1), const.cp(range2,2)];
        cNeg    = [const.cn(range1,1), const.cn(range2,2)];
        
        switch branch
            %Main branches:
            case 'p11'
                pref 	= (J2-lambda-0.5).*(J2+lambda+0.5) ./ (4*J2.*cNeg(:,1).*cNeg(:,2));
                out 	= pref.*( uNeg(:,1).*uNeg(:,2) + 4*(J2-lambda+0.5).*(J2+lambda-0.5) ).^2;
                
            case 'p22'
                pref 	= (J2-lambda-0.5).*(J2+lambda+0.5) ./ (4*J2.*cPos(:,1).*cPos(:,2));
                out 	= pref.*( uPos(:,1).*uPos(:,2) + 4*(J2-lambda+0.5).*(J2+lambda-0.5) ).^2;
                
            case 'q11'
                pref 	= (J2+0.5) ./ (2*J2.*(J2+1).*cNeg(:,1).*cNeg(:,2));
                out 	= pref.*( (lambda+0.5).*uNeg(:,1).*uNeg(:,2) ...
                    +  	4*(lambda-0.5).*(J2-lambda+0.5).*(J2+lambda+0.5) ).^2;
                
            case 'q22'
                pref 	= (J2+0.5) ./ (2*J2.*(J2+1).*cPos(:,1).*cPos(:,2));
                out 	= pref.*( (lambda+0.5).*uPos(:,1).*uPos(:,2) ...
                    +  	4*(lambda-0.5).*(J2-lambda+0.5).*(J2+lambda+0.5) ).^2;
                
            case 'r11'
                pref 	= (J1-lambda-0.5).*(J1+lambda+0.5) ./ (4*J1.*cNeg(:,1).*cNeg(:,2));
                out 	= pref.*( uNeg(:,1).*uNeg(:,2) + 4*(J1-lambda+0.5).*(J1+lambda-0.5) ).^2;
                
            case 'r22'
                pref 	= (J1-lambda-0.5).*(J1+lambda+0.5) ./ (4*J1.*cPos(:,1).*cPos(:,2));
                out 	= pref.*( uPos(:,1).*uPos(:,2) + 4*(J1-lambda+0.5).*(J1+lambda-0.5) ).^2;
                
                
                %satallite branches:
            case 'p21'
                pref 	= (J2-lambda-0.5).*(J2+lambda+0.5) ./ (4*J2.*cPos(:,1).*cNeg(:,2));
                out 	= pref.*( uPos(:,1).*uNeg(:,2) - 4*(J2-lambda+0.5).*(J2+lambda-0.5) ).^2;
                
            case 'p12'
                pref 	= (J2-lambda-0.5).*(J2+lambda+0.5) ./ (4*J2.*cNeg(:,1).*cPos(:,2));
                out 	= pref.*( uNeg(:,1).*uPos(:,2) - 4*(J2-lambda+0.5).*(J2+lambda-0.5) ).^2;
                
            case 'q21'
                pref 	= (J2+0.5) ./ (2*J2.*(J2+1).*cPos(:,1).*cNeg(:,2));
                out 	= pref.*( (lambda+0.5).*uPos(:,1).*uNeg(:,2) ...
                    -  	4*(lambda-0.5).*(J2-lambda+0.5).*(J2+lambda+0.5) ).^2;
                
            case 'q12'
                pref 	= (J2+0.5) ./ (2*J2.*(J2+1).*cNeg(:,1).*cPos(:,2));
                out 	= pref.*( (lambda+0.5).*uNeg(:,1).*uPos(:,2) ...
                    -  	4*(lambda-0.5).*(J2-lambda+0.5).*(J2+lambda+0.5) ).^2;
                
            case 'r21'
                pref 	= (J1-lambda-0.5).*(J1+lambda+0.5) ./ (4*J1.*cPos(:,1).*cNeg(:,2));
                out 	= pref.*( uPos(:,1).*uNeg(:,2) - 4*(J1-lambda+0.5).*(J1+lambda-0.5) ).^2;
                
            case 'r12'
                pref 	= (J1-lambda-0.5).*(J1+lambda+0.5) ./ (4*J1.*cNeg(:,1).*cPos(:,2));
                out 	= pref.*( uNeg(:,1).*uPos(:,2) - 4*(J1-lambda+0.5).*(J1+lambda-0.5) ).^2;
                
        end
    end

    function [out] = hlTripletCalc(branch, const, lambda, Y1, Y2, J1, J2, range1, range2)
        J1 		= J1(range1,1);
        J2 		= J2(range2,1);
        u1Pos 	= [const.u1p(range1,1), const.u1p(range2,2)];
        u1Neg 	= [const.u1n(range1,1), const.u1n(range2,2)];
        u3Pos 	= [const.u3p(range1,1), const.u3p(range2,2)];
        u3Neg 	= [const.u3n(range1,1), const.u3n(range2,2)];
        c1 		= [const.c1(range1,1), const.c1(range2,2)];
        c2 		= [const.c2(range1,1), const.c2(range2,2)];
        c3 		= [const.c3(range1,1), const.c3(range2,2)];
        
        
        switch branch
            %main branches
            case 'p11'
                pref    = 1*(J2-lambda).*(J2+lambda) .* 1./(16*J2.*c1(:,1).*c1(:,2));
                out 	= pref.*( (J2-lambda+1).*(J2+lambda-1) .* u1Pos(:,1).*u1Pos(:,2) ...
                            +     (J2-lambda-1).*(J2+lambda+1) .* u1Neg(:,1).*u1Neg(:,2) ...
                            + 	8*(J2-lambda-1).*(J2+lambda-1) .* (J2-lambda).*(J2+lambda) ).^2;
                
            case 'p22'
                pref    = 4*(J2-lambda).*(J2+lambda) .* 1./(J2.*c2(:,1).*c2(:,2));
                out 	= pref.*( (0.5*lambda^2*(Y1-2)*(Y2-2)) ...
                            + 	(J2-lambda-1).*(J2+lambda+1) ...
                            + 	(J2-lambda+1).*(J2+lambda-1) ).^2;
                
            case 'p33'
                pref    = 1*(J2-lambda).*(J2+lambda) .* 1./(16*J2.*c3(:,1).*c3(:,2));
                out 	= pref.*( (J2-lambda+1).*(J2+lambda-1) .* u3Neg(:,1).*u3Neg(:,2) ...
                            +	  (J2-lambda-1).*(J2+lambda+1) .* u3Pos(:,1).*u3Pos(:,2) ...
                            +	8*(J2-lambda+1).*(J2+lambda+1) .*(J2-lambda).*(J2+lambda) ).^2;
                
            case 'q11'
                pref    = 1*(2*J2+1) .* 1./(16*J2.*(J2+1).*c1(:,1).*c1(:,2));
                out  	= pref.*( (lambda-1).*(J2-lambda+1).*(J2+lambda)  .* u1Pos(:,1).*u1Pos(:,2) ...
                            +     (lambda+1).*(J2-lambda)  .*(J2+lambda+1).* u1Neg(:,1).*u1Neg(:,2) ...
                            + 	 8*lambda*(J2-lambda).^2 .*(J2+lambda).^2 ).^2;
                
            case 'q22'
                pref    = 4*(2*J2+1) .* 1./(J2.*(J2+1).*c2(:,1).*c2(:,2));
                out  	= pref.*( (0.5* lambda^3 * (Y1-2)*(Y2-2)) ...
                            + 	(lambda+1).*(J2-lambda)  .*(J2+lambda+1) ...
                            + 	(lambda-1).*(J2-lambda+1).*(J2+lambda) ).^2;
                
            case 'q33'
                pref    = 1*(2*J2+1) .* 1./(16*J2.*(J2+1).*c3(:,1).*c3(:,2));
                out  	= pref.*( (lambda-1).*(J2-lambda+1).*(J2+lambda)  .* u3Neg(:,1).*u3Neg(:,2) ...
                            +     (lambda+1).*(J2-lambda)  .*(J2+lambda+1).* u3Pos(:,1).*u3Pos(:,2) ...
                            + 	 8*lambda*(J2-lambda+1).^2 .*(J2+lambda+1).^2 ).^2;
                
            case 'r11'
                pref    = 1*(J1-lambda).*(J1+lambda) .* 1./(16*J1.*c1(:,1).*c1(:,2));
                out   	= pref.*( (J1-lambda+1).*(J1+lambda-1) .* u1Pos(:,1).*u1Pos(:,2) ...
                            +     (J1-lambda-1).*(J1+lambda+1) .* u1Neg(:,1).*u1Neg(:,2) ...
                            + 	8*(J2-lambda+1).*(J2+lambda+1).*(J2-lambda).*(J2+lambda) ).^2;
                
            case 'r22'
                pref 	= 4*(J1-lambda).*(J1+lambda) .* 1./(J1.*c2(:,1).*c2(:,2));
                out 	= pref.*( (0.5*lambda^2*(Y1-2)*(Y2-2)) ...
                            + 	  (J1-lambda-1).*(J1+lambda+1) ...
                            + 	  (J1-lambda+1).*(J1+lambda-1) ).^2;
                
            case 'r33'
                pref  	= 1*(J1-lambda).*(J1+lambda) .* 1./(16*J1.*c3(:,1).*c3(:,2));
                out 	= pref.*( (J1-lambda+1).*(J1+lambda-1).* u3Neg(:,1).*u3Neg(:,2) ...
                            +     (J1-lambda-1).*(J1+lambda+1).* u3Pos(:,1).*u3Pos(:,2) ...
                            +   8*(J2-lambda+1).*(J2-lambda+2).*(J2+lambda+1).*(J2+lambda+ 2) ).^2;
                
            %satallite branches
            %p-branches
            case 'p12'
                pref    = 1*(J2-lambda).*(J2+lambda) .* 1./(2*J2.*c1(:,1).*c2(:,2));
                out 	= pref.*( (J2-lambda+1).*(J2+lambda-1) .* u1Pos(:,1) ...
                            -     (J2-lambda-1).*(J2+lambda+1) .* u1Neg(:,1) ...
                            - 	2*lambda.*(Y2-2).*(J1-lambda).*(J1+lambda) ).^2;
                
            case 'p13'
                pref    = 1*(J2-lambda).*(J2+lambda) .* 1./(16*J2.*c1(:,1).*c3(:,2));
                out 	= pref.*( (J2-lambda+1).*(J2+lambda-1) .* u1Pos(:,1).*u3Neg(:,2) ...
                            +     (J2-lambda-1).*(J2+lambda+1) .* u1Neg(:,1).*u3Pos(:,2) ...
                            - 	8*(J1-lambda).*(J2-lambda+1).*(J1+lambda).*(J2+lambda+1) ).^2;
                
            case 'p21'
                pref    = 1*(J2-lambda).*(J2+lambda) .* 1./(2*J2.*c2(:,1).*c1(:,2));
                out 	= pref.*( (J2-lambda+1).*(J2+lambda-1) .* u1Pos(:,2) ...
                            -     (J2-lambda-1).*(J2+lambda+1) .* u1Neg(:,2) ...
                            - 	2*lambda.*(Y1-2).*(J2-lambda).*(J2+lambda) ).^2;
                
            case 'p23'
                pref    = 1*(J2-lambda).*(J2+lambda) .* 1./(2*J2.*c2(:,1).*c3(:,2));
                out 	= pref.*( (J2-lambda+1).*(J2+lambda-1) .* u3Neg(:,2) ...
                            -     (J2-lambda-1).*(J2+lambda+1) .* u3Pos(:,2) ...
                            + 	2*lambda.*(Y1-2).*(J2-lambda+1).*(J2+lambda+1) ).^2;
                
            case 'p31'
                pref    = 1*(J2-lambda).*(J2+lambda) .* 1./(16*J2.*c3(:,1).*c1(:,2));
                out 	= pref.*( (J2-lambda+1).*(J2+lambda-1) .* u3Neg(:,1).*u1Pos(:,2) ...
                            +     (J2-lambda-1).*(J2+lambda+1) .* u3Pos(:,1).*u1Neg(:,2) ...
                            - 	8*(J2-lambda).^2.*(J2+lambda).^2 ).^2;
                
            case 'p32'
                pref    = 1*(J2-lambda).*(J2+lambda) .* 1./(2*J2.*c3(:,1).*c2(:,2));
                out 	= pref.*( (J2-lambda+1).*(J2+lambda-1) .* u3Neg(:,1) ...
                            -     (J2-lambda-1).*(J2+lambda+1) .* u3Pos(:,1) ...
                            + 	2*lambda.*(Y2-2).*(J2-lambda).*(J2+lambda) ).^2;
            %q-branches
            case 'q12'
                pref    = 1*(2*J2+1) .* 1./(2*J2.*(J2+1).*c1(:,1).*c2(:,2));
                out  	= pref.*( (lambda-1).*(J2-lambda+1).*(J2+lambda)  .* u1Pos(:,1) ...
                            -     (lambda+1).*(J2-lambda)  .*(J2+lambda+1).* u1Neg(:,1) ...
                            - 	 2*lambda.^2*(J2-lambda).*(J2+lambda).*(Y2-2) ).^2;
                
            case 'q13'
                pref    = 1*(2*J2+1) .* 1./(16*J2.*(J2+1).*c1(:,1).*c3(:,2));
                out  	= pref.*( (lambda-1).*(J2-lambda+1).*(J2+lambda)  .* u1Pos(:,1).*u3Neg(:,2) ...
                            +     (lambda+1).*(J2-lambda)  .*(J2+lambda+1).* u1Neg(:,1).*u3Pos(:,2) ...
                            - 	 8*lambda*(J2-lambda).*(J2-lambda+1).*(J2+lambda).*(J2+lambda+1) ).^2;
                
            case 'q21'
                pref    = 1*(2*J2+1) .* 1./(2*J2.*(J2+1).*c2(:,1).*c1(:,2));
                out  	= pref.*( (lambda-1).*(J2-lambda+1).*(J2+lambda)  .* u1Pos(:,2) ...
                              -     (lambda+1).*(J2-lambda)  .*(J2+lambda+1).* u1Neg(:,2) ...
                              - 	 2*lambda^2*(J2-lambda).*(J2+lambda).*(Y1-2)).^2;
                
            case 'q23'
                pref    = 1*(2*J2+1) .* 1./(2*J2.*(J2+1).*c2(:,1).*c3(:,2));
                out  	= pref.*( (lambda-1).*(J2-lambda+1).*(J2+lambda)  .* u3Neg(:,2) ...
                            -     (lambda+1).*(J2-lambda)  .*(J2+lambda+1).* u3Pos(:,2) ...
                            + 	 2*lambda^2*(J2-lambda+1).*(J2+lambda+1).*(Y1-2) ).^2;
                
            case 'q31'
                pref    = 1*(2*J2+1) .* 1./(16*J2.*(J2+1).*c3(:,2).*c1(:,2));
                out  	= pref.*( (lambda-1).*(J2-lambda+1).*(J2+lambda)  .* u3Neg(:,1).*u1Pos(:,2) ...
                            +     (lambda+1).*(J2-lambda)  .*(J2+lambda+1).* u3Pos(:,1).*u1Neg(:,2) ...
                            - 	 8*lambda*(J2-lambda).*(J2-lambda+1).*(J2+lambda).*(J2+lambda+1) ).^2;
                
            case 'q32'
                pref    = 1*(2*J2+1) .* 1./(2*J2.*(J2+1).*c3(:,1).*c2(:,2));
                out  	= pref.*( (lambda-1).*(J2-lambda+1).*(J2+lambda)  .* u3Neg(:,1) ...
                            -     (lambda+1).*(J2-lambda)  .*(J2+lambda+1).* u3Pos(:,1) ...
                            + 	 2*lambda^2*(J2-lambda+1).*(J2+lambda+1).*(Y2-2) ).^2;
                
            %R-branch
            case 'r12'
                pref    = 1*(J1-lambda).*(J1+lambda) .* 1./(2*J1.*c1(:,1).*c2(:,2));
                out 	= pref.*( (J1-lambda+1).*(J1+lambda-1) .* u1Pos(:,1) ...
                            -     (J1-lambda-1).*(J1+lambda+1) .* u1Neg(:,1) ...
                            - 	2*lambda.*(Y2-2).*(J1-lambda).*(J1+lambda) ).^2;
                
            case 'r13'
                pref    = 1*(J1-lambda).*(J1+lambda) .* 1./(16*J1.*c1(:,1).*c3(:,2));
                out 	= pref.*( (J1-lambda+1).*(J1+lambda-1) .* u1Pos(:,1).*u3Neg(:,2) ...
                            +     (J1-lambda-1).*(J1+lambda+1) .* u1Neg(:,1).*u3Pos(:,2) ...
                            - 	8*(J1-lambda).^2 .* (J1+lambda).^2 ).^2;
                
            case 'r21'
                pref    = 1*(J1-lambda).*(J1+lambda) .* 1./(2*J1.*c2(:,1).*c1(:,2));
                out 	= pref.*( (J1-lambda+1).*(J1+lambda-1) .* u1Pos(:,2) ...
                            -     (J1-lambda-1).*(J1+lambda+1) .* u1Neg(:,2) ...
                            - 	2*lambda.*(Y1-2).*(J2-lambda).*(J2+lambda) ).^2;
                
            case 'r23'
                pref    = 1*(J1-lambda).*(J1+lambda) .* 1./(2*J1.*c2(:,1).*c3(:,2));
                out 	= pref.*( (J1-lambda+1).*(J1+lambda-1) .* u3Neg(:,2) ...
                            -     (J1-lambda-1).*(J1+lambda+1) .* u3Pos(:,2) ...
                            + 	2*lambda.*(Y1-2).*(J2-lambda+1).*(J2+lambda+1) ).^2;
                
            case 'r31'
                pref    = 1*(J1-lambda).*(J1+lambda) .* 1./(16*J1.*c3(:,1).*c1(:,2));
                out 	= pref.*( (J1-lambda+1).*(J1+lambda-1) .* u3Neg(:,1).*u1Pos(:,2) ...
                            +     (J1-lambda-1).*(J1+lambda+1) .* u3Pos(:,1).*u1Neg(:,2) ...
                            - 	8*(J1-lambda-1).*(J1-lambda+1).*(J1+lambda-1).*(J1+lambda+1) ).^2;
                
            case 'r32'
                pref    = 1*(J1-lambda).*(J1+lambda) .* 1./(2*J1.*c3(:,1).*c2(:,2));
                out 	= pref.*( (J1-lambda+1).*(J1+lambda-1) .* u3Neg(:,1) ...
                            -     (J1-lambda-1).*(J1+lambda+1) .* u3Pos(:,1) ...
                            + 	2*lambda.*(Y2-2).*(J1-lambda+1).*(J1+lambda+1) ).^2;
                
        end
    end

    function [out] = hlTripletCalcDel1(branch, const, lambda, Y1, Y2, J1, J2, range1, range2)
        J1 		= J1(range1,1);
        J2 		= J2(range2,1);
        u1Pos 	= [const.u1p(range1,1), const.u1p(range2,2)];
        u1Neg 	= [const.u1n(range1,1), const.u1n(range2,2)];
        u3Pos 	= [const.u3p(range1,1), const.u3p(range2,2)];
        u3Neg 	= [const.u3n(range1,1), const.u3n(range2,2)];
        c1 		= [const.c1(range1,1), const.c1(range2,2)];
        c2 		= [const.c2(range1,1), const.c2(range2,2)];
        c3 		= [const.c3(range1,1), const.c3(range2,2)];
        
        switch branch
            %main branches
            case 'p11'
		        pref 	= (J1-lambda).*(J2-lambda)./(32*J2.*c1(:,1).c1(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J1+lambda+1) .* u1Pos(:,1).*u1Pos(:,2) ...
					        + 	  (J2-lambda-2).*(J1+lambda+2) .* u1Neg(:,1).*u1Neg(:,2) ...
					        + 	8*(J2-lambda-2).*(J2-lambda).*(J2+lambda).^2 ).^2;
		        
                
            case 'p22'
		        pref 	= 2*(J1-lambda).*(J2-lambda)./(J2.*c2(:,1).*c2(:,2));
		        out 	= pref.*(0.5*lambda*(lambda+1)*(Y1-2)*(Y2-2) ...
					        + 	  (J2-lambda+1).*(J1+lambda+1) ...
					        + 	  (J2-lambda-2).*(J1+lambda+2) ).^2;
                
            case 'p33'
                pref 	= (J1-lambda).*(J2-lambda)./(32*J2.*c3(:,1).c3(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J1+lambda+1) .* u3Neg(:,1).*u3Neg(:,2) ...
					        + 	  (J2-lambda-2).*(J1+lambda+2) .* u3Pos(:,1).*u3Pos(:,2) ...
					        + 	8*(J2-lambda-1).*(J2-lambda+1).*(J2+lambda+1).^2 ).^2;
		        
                
            case 'q11'
		        pref 	= (J2-lambda).*(J2+lambda+1).*(2*J2+1)./(32*J2.*(J2+1).*c1(:,1).*c1(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J2+lambda)   .* u1Pos(:,1).*u1Pos(:,2) ...
					        +	  (J2-lambda-1).*(J2+lambda+2) .* u1Neg(:,1).*u1Neg(:,2) ...
					        + 	8*(J2-lambda-1).*(J2-lambda).*(J2+lambda).*(J2+lambda+1) ).^2;
                
            case 'q22'
		        pref 	= 2*(J2-lambda).*(J2+lambda+1).*(2*J2+1)./(J2.*(J2+1).*c2(:,1).*c2(:,2));
		        out 	= pref.*(0.5*lambda*(lambda+1)*(Y1-2)*(Y2-2) ...
					        + 	  (J2-lambda+1).*(J2+lambda) ...
					        + 	  (J2-lambda-1).*(J2+lambda+2) ).^2;
        
        
            case 'q33'
			        pref 	= (J2-lambda).*(J2+lambda+1).*(2*J2+1)./(32*J2.*(J2+1).*c3(:,1).*c3(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J2+lambda)   .* u3Neg(:,1).*u3Neg(:,2) ...
					        +	  (J2-lambda-1).*(J2+lambda+2) .* u3Pos(:,1).*u3Pos(:,2) ...
					        + 	8*(J2-lambda).*(J2-lambda+1).*(J2+lambda+1).*(J2+lambda+2) ).^2;
                
                
            case 'r11'
		        pref 	= (J1+lambda).*(J2+lambda+2)./(32*J1.*c1(:,1).c1(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J1+lambda-1) .* u1Pos(:,1).*u1Pos(:,2) ...
					        + 	  (J2-lambda)  .*(J1+lambda+2) .* u1Neg(:,1).*u1Neg(:,2) ...
					        + 	8*(J2+lambda+2).*(J2+lambda).*(J2-lambda).^2 ).^2;
                
            case 'r22'
		        pref 	= 2*(J1+lambda).*(J2+lambda+2)./(J1.*c2(:,1).*c2(:,2));
		        out 	= pref.*(0.5*lambda*(lambda+1)*(Y1-2)*(Y2-2) ...
					        + 	  (J2-lambda+1).*(J1+lambda-1) ...
					        + 	  (J2-lambda)  .*(J1+lambda+2) ).^2;
                
            case 'r33'
		        pref 	= (J1+lambda).*(J2+lambda+2)./(32*J1.*c3(:,1).c3(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J1+lambda-1) .* u3Neg(:,1).*u3Neg(:,2) ...
					        + 	  (J2-lambda)  .*(J1+lambda+2) .* u3Pos(:,1).*u3Pos(:,2) ...
					        + 	8*(J2+lambda+3).*(J2+lambda+1).*(J2-lambda+1).^2 ).^2;
                
            %satallite branches
            %p-branches
            case 'p12'
		        pref 	= (J1-lambda).*(J2-lambda)./(4*J2.*c1(:,1).c2(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J1+lambda+1) .* u1Pos(:,1) ...
					        - 	  (J2-lambda-2).*(J1+lambda+2) .* u1Neg(:,1) ...
					        - 	2*lambda*(J2-lambda-2).*(J2+lambda)*(Y2-2) ).^2;
					        
                
            case 'p13'
		        pref 	= (J1-lambda).*(J2-lambda)./(32*J2.*c1(:,1).c3(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J1+lambda+1) .* u1Pos(:,1).*u3Neg(:,2) ...
					        + 	  (J2-lambda-2).*(J1+lambda+2) .* u1Neg(:,1).*u3Pos(:,2) ...
					        - 	8*(J2-lambda-2).*(J2-lambda+1).*(J2+lambda).*(J2+lambda+1) ).^2;
                
            case 'p21'
		        pref 	= (J1-lambda).*(J2-lambda)./(4*J2.*c2(:,1).c1(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J1+lambda+1) .* u1Pos(:,2) ...
					        - 	  (J2-lambda-2).*(J1+lambda+2) .* u1Neg(:,2) ...
					        - 	2*(lambda+1)*(J2-lambda).*(J2+lambda)*(Y1-2) ).^2;
                
            case 'p23'
		        pref 	= (J1-lambda).*(J2-lambda)./(4*J2.*c2(:,1).c3(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J1+lambda+1) .* u3Neg(:,2) ...
					        - 	  (J2-lambda-2).*(J1+lambda+2) .* u3Pos(:,2) ...
					        + 	2*(lambda+1)*(J2-lambda+1).*(J2+lambda+1)*(Y1-2) ).^2;
                
            case 'p31'
		        pref 	= (J1-lambda).*(J2-lambda)./(32*J2.*c3(:,1).c1(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J1+lambda+1) .* u3Neg(:,1).*u1Pos(:,2) ...
					        + 	  (J2-lambda-2).*(J1+lambda+2) .* u3Pos(:,1).*u1Neg(:,2) ...
					        - 	8*(J2-lambda-1).*(J2-lambda).*(J2+lambda).*(J2+lambda+1) ).^2;
                
            case 'p32'
		        pref 	= (J1-lambda).*(J2-lambda)./(4*J2.*c3(:,1).c2(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J1+lambda+1) .* u3Neg(:,1) ...
					        - 	  (J2-lambda-2).*(J1+lambda+2) .* u3Pos(:,1) ...
					        + 	2*lambda*(J2-lambda-1).*(J2+lambda+1)*(Y2-2) ).^2;
            %q-branches
            case 'q12'
		        pref 	= (J2-lambda).*(J2+lambda+1).*(2*J2+1)./(4*J2.*(J2+1).*c1(:,1).*c2(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J2+lambda)   .* u1Pos(:,1) ...
					        -	  (J2-lambda-1).*(J2+lambda+2) .* u1Neg(:,1) ...
					        - 	2*lambda*(J2-lambda-1).*(J2+lambda+1).*(Y2-2) ).^2;
                
            case 'q13'
		        pref 	= (J2-lambda).*(J2+lambda+1).*(2*J2+1)./(32*J2.*(J2+1).*c1(:,1).*c3(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J2+lambda)   .* u1Pos(:,1).*u3Neg(:,2) ...
					        +	  (J2-lambda-1).*(J2+lambda+2) .* u1Neg(:,1).*u3Pos(:,2) ...
					        - 	8*(J2-lambda-1).*(J2-lambda+1).*(J2+lambda+1).^2 ).^2;
                
            case 'q21'
		        pref 	= (J2-lambda).*(J2+lambda+1).*(2*J2+1)./(4*J2.*(J2+1).*c2(:,1).*c1(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J2+lambda)   .* u1Pos(:,2) ...
					        -	  (J2-lambda-1).*(J2+lambda+2) .* u1Neg(:,2) ...
					        - 	2*(lambda+1)*(J2-lambda).*(J2+lambda).*(Y1-2) ).^2;
                
            case 'q23'
		        pref 	= (J2-lambda).*(J2+lambda+1).*(2*J2+1)./(4*J2.*(J2+1).*c2(:,1).*c3(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J2+lambda)   .* u3Neg(:,2) ...
					        -	  (J2-lambda-1).*(J2+lambda+2) .* u3Pos(:,2) ...
					        + 	2*(lambda+1)*(J2-lambda+1).*(J2+lambda+1).*(Y1-2) ).^2;
                
            case 'q31'
		        pref 	= (J2-lambda).*(J2+lambda+1).*(2*J2+1)./(32*J2.*(J2+1).*c3(:,1).*c1(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J2+lambda)   .* u3Neg(:,1).*u1Pos(:,2) ...
					        +	  (J2-lambda-1).*(J2+lambda+2) .* u3Pos(:,1).*u1Neg(:,2) ...
					        - 	8*(J2+lambda+2).*(J2+lambda).*(J2-lambda).^2 ).^2;
                
            case 'q32'
		        pref 	= (J2-lambda).*(J2+lambda+1).*(2*J2+1)./(4*J2.*(J2+1).*c3(:,1).*c2(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J2+lambda)   .* u3Neg(:,1) ...
					        -	  (J2-lambda-1).*(J2+lambda+2) .* u3Pos(:,1) ...
					        + 	2*lambda*(J2-lambda).*(J2+lambda+2).*(Y2-2) ).^2;
                
            %R-branch
            case 'r12'
		        pref 	= (J1+lambda).*(J2+lambda+2)./(4*J1.*c1(:,1).c2(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J1+lambda-1) .* u1Pos(:,1) ...
					        - 	  (J2-lambda)  .*(J1+lambda+2) .* u1Neg(:,1) ...
					        - 	2*lambda*(J2-lambda).*(J2+lambda+2).*(Y2-2) ).^2;
                
            case 'r13'
		        pref 	= (J1+lambda).*(J2+lambda+2)./(32*J1.*c1(:,1).c3(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J1+lambda-1) .* u1Pos(:,1).*u3Neg(:,2) ...
					        + 	  (J2-lambda)  .*(J1+lambda+2) .* u1Neg(:,1).*u3Pos(:,2) ...
					        - 	8*(J2-lambda).*(J2-lambda+1).*(J2+lambda+1).*(J2+lambda+2) ).^2;
                
            case 'r21'
		        pref 	= (J1+lambda).*(J2+lambda+2)./(4*J1.*c2(:,1).c1(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J1+lambda-1) .* u1Pos(:,2) ...
					        - 	  (J2-lambda)  .*(J1+lambda+2) .* u1Neg(:,2) ...
					        - 	2*(lambda+1)*(J2-lambda).*(J2+lambda).*(Y1-2) ).^2;
                
            case 'r23'
		        pref 	= (J1+lambda).*(J2+lambda+2)./(4*J1.*c2(:,1).c3(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J1+lambda-1) .* u3Neg(:,2) ...
					        - 	  (J2-lambda)  .*(J1+lambda+2) .* u3Pos(:,2) ...
					        + 	2*(lambda+1)*(J2-lambda+1).*(J2+lambda+1).*(Y1-2) ).^2;
                
            case 'r31'
		        pref 	= (J1+lambda).*(J2+lambda+2)./(32*J1.*c3(:,1).c1(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J1+lambda-1) .* u3Neg(:,1).*u1Pos(:,2) ...
					        + 	  (J2-lambda)  .*(J1+lambda+2) .* u3Pos(:,1).*u1Neg(:,2) ...
					        - 	8*(J2-lambda).*(J2-lambda+1).*(J2+lambda).*(J2+lambda+3) ).^2;
                
            case 'r32'
		        pref 	= (J1+lambda).*(J2+lambda+2)./(4*J1.*c3(:,1).c2(:,2));
		        out 	= pref.*( (J2-lambda+1).*(J1+lambda-1) .* u3Neg(:,1) ...
					        - 	  (J2-lambda)  .*(J1+lambda+2) .* u3Pos(:,1) ...
					        + 	2*lambda*(J2-lambda+1).*(J2+lambda+3).*(Y2-2) ).^2;
        end
	end

end