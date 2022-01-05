function [] = OptimizeMolParameters()
    %syst consnants (user input)
    syst  = '1-';
    V1      = [0:4,1:2];
    V2      = [0:4,0,1];
    Jf      = 80;
    filename    = 'FNS-4000Tr-4000Tv-Specair_0.08.txt';
	
    %Spectroscopic paramters (fix)
    spec.fwhm = 0.08;           %fwhm of the spectrometer [nm]
    spec.res  = spec.fwhm/10;   %spectrometer resolution (wavelenght/bin) [nm]
    spec.p    = 1;              %relative magnitude of gaussian function.
    spec.wl   = (355:spec.res:391.6)';    %wavelength range of interest.
	
	%Temperature (user input)
    temp.rot    = 4000;
    temp.vib    = 4000;

    %measured data (user input)
    measData    = importdata(filename);
    measWl      = measData(:,1);
    measI       = measData(:,2);

    %Define Initial points and upper/lower bounds
    x0  = ones(2,21);
    [ub,lb]    = getConstraints(x0);

    y0  = ones(2,length(V1));
    y0(2,:)     = 0;

    lby = y0*0.8; 
    uby = y0*2; 
    lby(:,1) = y0(:,1);
    uby(:,1) = y0(:,1);

    x0  = [x0, y0];
    lb  = [lb, lby];
    ub  = [ub, uby];
	
	%Plot data with intial guess
    err0  = probFunction(x0, syst, V1, V2, Jf, spec, temp, measI, measWl,1); 
	display(err0);
    title('Before Corrections');

	%Start optimization fucntion
    ms  = MultiStart('UseParallel', true);
    options = optimoptions('fmincon','FunctionTolerance',1e-6);
    problem = createOptimProblem('fmincon','objective',@(x0) probFunction(x0, syst, V1, V2, Jf, spec, temp, measI, measWl,0),...
                                    'x0', x0, 'lb',lb,'ub',ub,'options',options);   % Create optimization problem to solve
    tic;
    out = run(ms,problem,4);
    display(toc)
    
	%Plot data with optimized paratmeters
    err1  = probFunction(out, syst, V1, V2, Jf, spec, temp, measI, measWl,1); 
	display(err1);
    title('After Corrections');
    
	%Display other outputs
    fc = getFCFs(syst, V1, V2);
    fcCorrection  = out(1,22:end); display(fcCorrection);
    
    dhCorrection = out(:,1:21);
    [~,~] =    molecularConstants(syst, max([V1, V2]), dhCorrection);    

    function [err] = probFunction (factors, syst, V1, V2, Jf, spec, temp, measI, measWl, plotFigure)
        duncorr = factors(:,1:21);
		fcCorr  = factors(1,22:end);

        I       = zeros(size(spec.wl));
        
		const 	= molecularConstants(syst, max([V1, V2]), duncorr);
		
		for i   = 1:length(V1)
			dataNew     = waveNumber(syst, V1(i), V2(i), Jf, const);
			[dataNew.Sp, dataNew.Sq, dataNew.Sr]  = honlLondon(dataNew);
			dataNew.fc  = frankCondon(dataNew)*fcCorr(i);
			[wl, inten]         = simulateSpectrum(spec, dataNew, temp);
			I                   = I + inten;
		end
        
        simWl   = wl;
        simI    = I/max(I);

        if plotFigure == 1
            figure; hold on;
            plot(measWl, measI,'k', 'DisplayName','Input Spectra','LineWidth',1.5);
            plot(simWl, simI, 'b', 'DisplayName','Simulated Spectra','LineWidth',1);
    
            xlabel('Wavelength [nm]');
            ylabel('Relative Intensity [a.u.]');
            grid on;
            legend('show')
        end

        measI   = interp1(measWl, measI, simWl);
        err     = sum(sqrt((measI - simI).^2));
    end

    function [out] = getFCFs(syst, V1, V2)
        FC   = frankCondon([], syst);
        out  = zeros(1,length(V1));
        for i = 1:length(V1)
                v1 = V1(i) + 1;
                v2 = V2(i) + 1;
                out(i) = FC(v1,v2);
        end
    end

    function [ub, lb] = getConstraints(x0)
        %con1/con2 variable allocation is as follows: 
        %1	2	3	4	5	6	7	8	9	10	11	12	13	14 (index)
        %T0	w	wxe	wye	wze	wae	B0	a	ax	ay	az	De	Y1	G1 (quantity)
        lb          = ones(size(x0));
        ub          = ones(size(x0));
		r           = [1;1];
		%upper bounds
		ub(:,1) 	= 1.2;
		ub(:,2:6) 	= r*[1.8, 1.8, 1.7, 1.6, 1.6];
		ub(:,7:11) 	= r*[1.8, 1.8, 1.8, 1.8, 1.8];
		ub(:,12:14) = r*[1.1, 1, 1];
		ub(:,15:17) = r*[1.1, 1.1, 1];
		ub(:,18:20) = r*[1.2, 1, 1];
		ub(:,21) 	= 1.2;
		
		%lower bounds
		lb(:,1) 	= 0.8;
		lb(:,2:6) 	= r*[0.5, 0.5, 0.5, 0.5, 0.5];
		lb(:,7:11) 	= r*[0.5, 0.5, 0.5, 0.5, 0.5];
		lb(:,12:14) = r*[0.7, 1, 1];
		lb(:,15:17) = r*[0.7, 1, 1];
		lb(:,18:20) = r*[0.7, 1, 1];
		lb(:,21) 	= 1.2;
    end
    
end