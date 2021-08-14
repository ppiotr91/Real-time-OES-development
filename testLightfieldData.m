clc; clear; close all;

fname       = '33100-PWM700_10-10Ar-0.8N2.csv';
rawData     = importdata(fname);

WAVELENGTH  = unique(rawData(:,1));
INTENSITY   = reshape(rawData(:,2), length(WAVELENGTH), []);
INTENSITY   = sum(INTENSITY,2);

