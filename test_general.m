clc; clear; close all;
syst        = 'CNvio';
V1 			= 4;
V2 			= 4;
Jf 			= 170;

temp.rot    = 4000;
temp.vib    = 4000;

spec.fwhm = 0.08;           		%fwhm of the spectrometer [nm]
spec.res  = spec.fwhm/10;          	%spectrometer resolution (wavelenght/bin) [nm]
spec.p    = 1;              		%relative magnitude of gaussian function.
spec.wl   = (385:spec.res:389)';    %wavelength range of interest.

[data, sim]		= simulateBand(syst, V1, V2, Jf, spec, temp);
[sim] 			= simulateSystem(syst, V1, V2, Jf, spec, temp);

%PLOT OPTIONS:
plotData(data,'fd');
plotData(data, 'hl');
plotData(data, 'sim', sim);

	function [data, sim] = simulateBand(syst, V1, V2, Jf, spec, temp, const)
		if nargin < 7
			const = molecularConstants(syst, max([V1, V2]) );
		end
		data                = waveNumber(syst,V1,V2,Jf, const);
		[data.Sp, data.Sq, data.Sr]  = honlLondon(data);
		data.fc             = frankCondon(data);
		[sim.wl, sim.i]     = simulateSpectrum(spec, data, temp);
	end

	function [sim] = simulateSystem(syst, v1, v2, Jf, spec, temp) 	
        v1      = 0:v1;
        v2      = 0:v2;
		[V1,V2] = meshgrid(v1,v2);
		V1      = reshape(V1,[],1);
		V2      = reshape(V2,[],1);
		
		I       = zeros(size(spec.wl));
		
		const 	= molecularConstants(syst, max([v1, v2]) );
		for i   = 1:length(V1)
			dataNew     = waveNumber(syst, V1(i), V2(i), Jf, const);
			[dataNew.Sp, dataNew.Sq, dataNew.Sr]  = honlLondon(dataNew);
			dataNew.fc  = frankCondon(dataNew);
			[wl, inten]         = simulateSpectrum(spec, dataNew, temp);
			I                   = I + inten;
		end
		
		sim.wl  	= wl;
		sim.i   	= I/max(I);
	end

	function [] 	= contourPlotting() %%%%%%%%%%%% needs work
		%CONTOUR1 plotting
		figure; hold on;
		d = generateContourData('contour1');
		contour(d.X, d.Y, d.Z.a2, '--', 'ShowText', 'on', 'LineWidth', 2, 'DisplayName', 'Area 2')
		contour(d.X, d.Y, d.Z.a1, 'ShowText', 'on', 'DisplayName', 'Area 1')
		legend ('show')
		grid on;
	end