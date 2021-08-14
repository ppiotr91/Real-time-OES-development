function [WAVELENGTH,INTENSITY] = reshapeLF(fname)
% Reshapes the raw data from Lightfield and outputs a single row of
% wavelenghts and intensities.
%   ADDTIONAL NOTES ON USAGE:
%   The raw data from LightField does not have pixels binned. Instead light
%   field will give the wavelength and the intensity of all the pixels in
%   the CCD array. This function bin the raw data to field all the unique
%   wavelength along with the summed intensities corresponding to these
%   wavelengths
%       INPUTS:-
%           - fName:    Filename of the .csv file containing the wavelength
%                       and intensity infromation.
%       OUTPUTS:-
%           - WAVELENGTH:   Wavelength in nm.
%           - INTENSITY:    Summed intensity.
rawData     = importdata(fname);

WAVELENGTH  = unique(rawData(:,1));
INTENSITY   = reshape(rawData(:,2), length(WAVELENGTH), []);
INTENSITY   = sum(INTENSITY,2);
end

