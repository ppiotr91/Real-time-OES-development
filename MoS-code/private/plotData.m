function [] = plotData(DATA, dataName, SPEC)
if nargin < 3
    SPEC = [];
end
figure; hold on;
switch dataName
    case 'fd'
        fortratDiagram(DATA.V1, DATA.V2, DATA.J2, DATA.wlP, DATA.wlQ, DATA.wlR, DATA.range);
    
    case 'hl'
        hlDiagram(DATA.Sp, DATA.Sq, DATA.Sr, DATA.J1, DATA.range);

    case 'fcMap'
        system  = DATA.syst;
        fc      = frankCondon([], system);
        fcColorMap(fc);

    case 'sim'
        plot(SPEC.wl, SPEC.i/max(SPEC.i));
        xlabel('Wavelenght [nm]')
        ylabel('Normalized Intensity')
        title('Simulated Spectrum')

    otherwise
        error('plotData: Invalid plot requested...');
end

    function [] = fortratDiagram(V1, V2, J2, wlP, wlQ, wlR, range)
        maxI    = size(wlP,2);
        maxJ    = size(wlP,2);
        %P-Branch:
        range1  = range.p.j1;
        range2  = range.p.j2;
        txtPre  = 'P_{%d%d}';
        SymbColor   = 'r';
        iterateOverBranch(J2, wlP, txtPre, SymbColor ,range1, range2, maxI, maxJ);
        
        range1  = range.q.j1;
        range2  = range.q.j2;
        txtPre  = 'Q_{%d%d}';
        SymbColor   = 'b';
        iterateOverBranch(J2, wlQ, txtPre, SymbColor ,range1, range2, maxI, maxJ);
        
        range1  = range.r.j1;
        range2  = range.r.j2;
        txtPre  = 'R_{%d%d}';
        SymbColor   = 'g';
        iterateOverBranch(J2, wlR, txtPre, SymbColor ,range1, range2, maxI, maxJ);
        
        xlabel('Wavelength (nm)');
        ylabel('Rotational quantum number, J''''');
        title(sprintf('Fortrat Diagram for (v'', v'''') = (%d,%d)',V1, V2));
        legend('show'); hold off;
        
            function [] = iterateOverBranch (J, wl, textPref, color, range1, range2, maxI, maxJ)
                symb    = {'d', '^', 'v'};
                for i   = 1:maxI
                    for j = 1:maxJ
                        if  sum(range1(:,j,i)) == 0; continue; end
                        if  sum(range2(:,j,i)) == 0; continue; end
                        r1  = range1(range1(:,j,i) ~= 0 , j,i);
                        r2  = range2(range2(:,j,i) ~= 0 , j,i);
                        if size(J) > 1
                            j2 = J(:,j); 
                        else
                            j2 = J; 
                        end
                        
                        txt     = sprintf(textPref,i,j);
                        scatter(wl(r1,j,i), j2(r2), color, 'filled', char(symb(i)), 'DisplayName', txt);
                    end
                end
            end
    end

    function [] = hlDiagram(Sp, Sq, Sr, J1, range)
        maxI = size(Sp,2);
        maxJ = size(Sp,2);
        
        range1  = range.p.j1;
        range2  = range.p.j2;
        txtPre     = 'P_{%d%d}'; linecolor = 'r';
        iterateOverBranch(J1, Sp, txtPre, linecolor,range1, range2, maxI, maxJ);
        
        range1  = range.q.j1;
        range2  = range.q.j2;
        txtPre     = 'Q_{%d%d}'; linecolor = 'b';
        iterateOverBranch(J1, Sq, txtPre, linecolor, range1, range2, maxI, maxJ);
        
        range1  = range.r.j1;
        range2  = range.r.j2;
        txtPre     = 'R_{%d%d}'; linecolor = 'g';
        iterateOverBranch(J1, Sr, txtPre, linecolor, range1, range2, maxI, maxJ);
        legend('show','location','northwest');
        xlabel('J"');
        ylabel('Honl-London factors');
        hold off;
        
            function [] = iterateOverBranch (J, S, textPref, color, range1, range2, maxI, maxJ)
                symb    = {'-','--','-.'};
                for i   = 1:maxI
                    for j = 1:maxJ
                        if  sum(range1(:,j,i)) == 0; continue; end
                        if  sum(range2(:,j,i)) == 0; continue; end
                        r1  = range1(range1(:,j,i) ~= 0 , j,i);
                        r2  = range2(range2(:,j,i) ~= 0 , j,i);
                        txt     = sprintf(textPref,i,j);
                        if size(J) > 1
                            j1 = J(:,j); 
                        else
                            j1 = J; 
                        end
                        notation = char(strcat(symb(i),color));
                        plot(j1(r1), S(r1,j,i), notation,'DisplayName', txt);
                    end
                end
            end
    end

    function [] = fcColorMap(Z)
        lenv1    = size(Z,1);
        lenv2    = size(Z,2);
        
        v1  = linspace(0,lenv1-1, lenv1);
        v2  = linspace(0,lenv2-1, lenv2);
        [Y,X]   = meshgrid(v2,v1);
        surface(X,Y,Z);
        colormap(winter)
        colorbar
        xlabel('Upper vibrational level');
        ylabel('Lower vibrational level');
    end

end

