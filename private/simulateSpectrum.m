function [wavelength, totInten] = simulateSpectrum(specParamters, DATA, temp)
rotT    = temp.rot;
vibT    = temp.vib;

%p-branch
S       = DATA.Sp;
wl      = DATA.wlP;
range1  = DATA.range.p.j1;
range2  = DATA.range.p.j2;
intenP  = iterateOverBranch(S,wl,DATA, rotT, vibT, specParamters, range1, range2);

%q-branch
S       = DATA.Sq;
wl      = DATA.wlQ;
range1  = DATA.range.q.j1;
range2  = DATA.range.q.j2;
intenQ  = iterateOverBranch(S,wl,DATA, rotT, vibT, specParamters, range1, range2);

%r-branch
S       = DATA.Sr;
wl      = DATA.wlR;
range1  = DATA.range.r.j1;
range2  = DATA.range.r.j2;
intenR  = iterateOverBranch(S,wl,DATA, rotT, vibT, specParamters, range1, range2);

totInten  = intenP+ intenR + intenQ;
totInten  = sum(totInten,2);
totInten  = sum(totInten,3);

wavelength  = specParamters.wl;

    function [sumIntensity]     = branchInten(wl, S, I, wavelength, spec)
        scale       = 1e10;
        I           = 1./wl.^4 .* S.*I*scale;
        intensity   = zeros(size(wavelength,1),size(wl,1));
        w           = spec.fwhm;
        gaussPre    = spec.p/w * sqrt(4 * log(2)/pi);
        lorentzPre  = (1-spec.p) * 2/pi * w;
        smallNum    = eps(1);
        
        for k   = 1:length(S)
            if wl(k) == 0 || I(k) == 0
                continue;
            end
            coef = -4 * log(2)/w^2 * (wavelength - wl(k)).^2;
            intensity(:,k)  = I(k) * gaussPre * (exp(coef)+smallNum) + ...
                                     lorentzPre * 1./(w^2 + 4 * (wavelength - wl(k)).^2);
        end

        if sum(intensity,2) == 0
            warning('simulateSpectrum/branchIntensity: Calculated intensity is zero');
        end

        sumIntensity     = sum(intensity,2);
    end

    function [inten]                 = iterateOverBranch(S,wl,DATA, rotT, vibT, specParamters, range1, range2)
        hck     = 1.986e-23/1.3807e-23;        %K/cm
        Gv      = DATA.Gv1;
        fcFac   = DATA.fc;
        specWL  = specParamters.wl;
        inten   = zeros(size(specWL,1),size(S,2),size(S,2));
        
        for m   = 1:size(S,2)
            for n = 1:size(S,2)
                if  sum(range1(:,n,m)) == 0; continue; end
                if  sum(range2(:,n,m)) == 0; continue; end
                r1  = range1(range1(:,n,m) ~= 0 , n,m);

                Fj  = DATA.F1(r1,m);
                I   = fcFac * exp(- hck * Fj/rotT) * exp(-hck*Gv/vibT);
                inten(:,n,m)   = branchInten(wl(r1,n,m), S(r1,n,m), I, specWL, specParamters);
            end
        end
    end

end