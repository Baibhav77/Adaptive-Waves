function [Tvis, SHGC, LSG] = heat(spectrum,emis)
%% solar function calculates heat-related characteristics of low-E glasses
%% INPUT
% spectrum - column-wise [lambda, T, R] (wavelength, transmittance, reflectance)
% emis - emissivity, a portion of absorbance converted to heat
%% OUTPUT
% Tvis - visible light transmittance
% SHGC - solar heat gain coefficient
% LSG - light-to-solar gain ratio
%% REFERENCE
% for details, see:
% McCluney, Ross (1996), 
% Fenestration Solar Gain Analysis, 
% Florida Solar Energy Center/University of Central Florida, 
% CiteSeerX 10.1.1.30.2472
% https://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=5FA45D124D679C8491E771C626052CE3?doi=10.1.1.30.2472&rep=rep1&type=pdf
%% LOADING DATA
load('AM15');  % solar illumination
load('Veye'); % human response eye function
%% Tvis
Tvis = interp1(spectrum(:,1), spectrum(:,2), Veye(:,1), 'spline');
AM15_vis = interp1(AM15(:,1), AM15(:,2), Veye(:,1), 'spline');
Tvis = sum(Tvis.*Veye(:,2).*AM15_vis)/sum(Veye(:,2).*AM15_vis); % Eq.(19)
%% SHGC (see Eqs.(12),(14))
lam_IR = linspace(300, 2000, 1700); % change the largest lambda to whatever available
AT = spectrum(:,2) + emis*(1-spectrum(:,2)-spectrum(:,3));
AM15_IR = interp1(AM15(:,1), AM15(:,2), lam_IR, 'spline');
AT_IR = interp1(spectrum(:,1), AT, lam_IR, 'spline');
SHGC = sum(AT_IR.*AM15_IR)/sum(AM15_IR);
%% LSG
LSG = Tvis/SHGC; % Eq.(20)
end