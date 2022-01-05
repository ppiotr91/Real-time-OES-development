function [contour] = generateContourData (savedDataType, spec)
%   RUN CONTOUR1: cont=generateContourData('contour1', spec);
if nargin < 2
    spec.fwhm = 0.1;            %fwhm of the spectrometer [nm]
    spec.res  = spec.fwhm/10;   %spectrometer resolution (wavelenght/bin) [nm]
    spec.p    = 1;              %relative magnitude of gaussian function.
    spec.wl   = (310:spec.res:340)';    %wavelength range of interest.
end

switch savedDataType
    case 'contour1'
        %User Inputs:
        TrotLower   = 2000;
        TrotUpper   = 8000;
        TvibLower   = 2000;
        TvibUpper   = 12000;
        nVib        = 20;
        nRot        = 25;
        system      = '2+';
        
        rotTemp     = linspace(TrotLower, TrotUpper, nRot)';
        vibTemp     = linspace(TvibLower, TvibUpper, nVib);
        [contROT, contVIB]  = meshgrid(rotTemp,vibTemp);
        
        contour.X       = contVIB;
        contour.Y       = contROT;
        area1           = zeros(size(contROT));
        area2           = zeros(size(contROT));
        
        %Function Inputs   
        cData.v1upper   = 4;
        cData.v2upper   = 4;
        bins            = [ 336.50  337.50
                            333.20  335.20
                            314.75  315.75]; %bins used: [0-0 bandHead; 0-0RotLines; 1-0BandHead]
        cData.bins      = bins;                
        cData.Tvib      = vibTemp;
        for m   = 1:length(rotTemp)
            cData.Trot      = rotTemp(m);
            cont1Data       = contour1(system,cData,spec);
            area            = cont1Data.area;
            
            area1(:,m)           = area(:,2)./area(:,1);
            area2(:,m)           = area(:,3)./area(:,1);
        end
        
        contour.Z.a1    = area1;
        contour.Z.a2    = area2;
        
    otherwise
        error('generateContorData: Invalid dataType request')
end


    function [OUT]      = contour1(system,cData,spec)
        %UserInputs:
        binBounds   = cData.bins;
        v1f         = cData.v1upper;
        v2f         = cData.v2upper;
        Tvib        = cData.Tvib;
        Trot        = cData.Trot;
        
        out     = zeros(length(Tvib), length(binBounds));

        for j = 1:length(Tvib)
            temp.rot    = Trot;
            temp.vib    = Tvib(j);
            [wl, I] = roVibronicSpectrum (system, v1f, v2f, spec, temp);
            output  = dissectSpectrum(wl,I, binBounds);
            out(j,:)    = output';
        end
        
        %Outputs:
        OUT.area    = out;
        OUT.bins    = binBounds;
        
    end

    function [wl, I]    = roVibronicSpectrum (system, V1f, V2f, spec, temp)
        v1  = 0:V1f;
        v2  = 0:V2f;
        [V1,V2] = meshgrid(v1,v2);
        V1      = reshape(V1,[],1);
        V2      = reshape(V2,[],1);
        
        I       = zeros(size(spec.wl));
        for i   = 1:length(V1)
            data                = waveNumber(system,V1(i),V2(i),120);
            [data.Sp, data.Sq, data.Sr]  = honlLondon(data);
            data.fc             = frankCondon(data);
            [wl, inten]         = simulateSpectrum(spec, data, temp);
            I                   = I + inten;
        end
    end

    function [binInten] = dissectSpectrum(wl, I, bins)
        binInten    = zeros(length(bins),1);
        for k   = 1:length(bins)
            binLower    = bins(k,1);
            binUpper    = bins(k,2);
            indexLower  = find(wl>binLower, 1 );
            indexUpper  = find(wl<binUpper, 1, 'last' );
            binInten(k) = sum(I(indexLower:indexUpper));
        end
    end

end

%fucntion template
%     function [] = FUNCTION NAME ()
% 
%     end