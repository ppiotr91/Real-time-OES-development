function [] = franckCondon2()
    phi1 = calcWaveFunctions ('N2-C3Pi');

    plot(phi1(4,:))

    function [out] = getConstants(state)
        switch state
            case 'N2-C3Pi'
                    we     = 2047.178;
                    wxe    = 28.445;
                    mu     = 7.00153720;
                    re     = 1.14869;
        end
        out.we    = we;
        out.wxe   = wxe;
        out.mu    = mu;
        out.re    = re;
    end

    function [PHI] = calcWaveFunctions (state)
        rmin    = 0.8;
        rmax    = 2;
        N       = 100;
        h       = 0.0039901772;     %amu*cm^2/s^2
            
        const   = getConstants(state);
        we      = const.we;
        wxe     = const.wxe;
        mu      = const.mu;
        re      = const.re;
        
        vmax    = round(0.5 * we/wxe);
        beta    = 0.2454 * (mu*wxe)^0.5;
        k       = 2*vmax;
        x       = linspace(rmin,rmax,N) - re;
        De      = we^2/(4*wxe);
        k1      = 4*pi*(2*mu*De)^0.5/(beta*h);
        d       = 0.5*k1;
        
        c1      = k*exp(-beta*x);
        PHI     = zeros(vmax+1, N);
        for v   = 0:30
            c2      = (k - 2*v - 1)*0.5;
            L       = laguerreL(v,c2,c1);
            phi     = (k*beta).^2 .* exp(-0.5*c1) .* c1.^c2 .* L;
            PHI(v+1,:)    = phi;
        end
        normal  = sum(PHI,2);
    end
    
end

