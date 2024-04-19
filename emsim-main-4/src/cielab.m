function [L, a, b] = cielab(lam,spectrum)
%% cielab function calculates L*a*b parameters for a given spectrum
%% INPUT
% lam - wavelength
% spectrum - spectrum of interest
% refl - true or false, whether it is a reflectance spectrum or not
%% OUTPUT
%- L* = coordinate of brightness, visible gamut range [0, 100] (L = 100
%       indicates diffuse white; specular white might be higher);
%- a* = coordinate of chromaticity for red-gree, gamut range [-128, 127];
%- b* = coordinate of chromaticity for yellow-blue, gamut range [-128, 127];
%% REFERENCES
% http://www.brucelindbloom.com/index.html?Equations.html
% https://scipython.com/blog/converting-a-spectrum-to-a-colour/
%% LOADING DATA
load('colMatching'); load('D65'); load("Veye");
%% REARRANGING INPUT
spectrum(lam<380) = [];
lam(lam<380) = [];
spectrum(lam>780) = [];
lam(lam>780) = [];
lamtab = transpose(380:1:780);
colMatching = interp1(lamtab,colMatching,lam);
D65 = interp1(lamtab,D65,lam);
%% CALCULATING L*a*b
% XYZn = [0.950489,1,1.08884]; % 2 deg observer
XYZn = [0.94811,1,1.07304]; % 10 deg observer
Par2 = D65.* colMatching(:,2);
ContrX = (spectrum.*D65).*colMatching(:,1);
ContrY = (spectrum.*D65).*colMatching(:,2);
ContrZ = (spectrum.*D65).*colMatching(:,3);
X = sum(ContrX)/sum(Par2);
Y = sum(ContrY)/sum(Par2);
Z = sum(ContrZ)/sum(Par2);
Ref_X = X/XYZn(1);
Ref_Y = Y/XYZn(2);
Ref_Z = Z/XYZn(3);
%
if Ref_X  > 0.008856 
    f_X = Ref_X^(1/3);
elseif Ref_X  <= 0.008856
    f_X = 7.787*Ref_X + 16/116;
end
%
if Ref_Y > 0.008856
    f_Y = Ref_Y^(1/3);
    L = 116*f_Y - 16;
elseif Ref_Y <= 0.008856
    L = 903.3*Ref_Y;
    f_Y = 7.787*Ref_Y + 16/116;
end
%
if Ref_Z  > 0.008856
    f_Z = Ref_Z^(1/3);
elseif Ref_Z <= 0.008856 
    f_Z = 7.787*Ref_Z + 16/116;
end
%
a = 500*(f_X - f_Y);
b = 200*(f_Y - f_Z);
%
%% RE-NORMALIZATION DUE TO VIRACON DATA
% if refl
%     L = 1.05*L;
%     a = 0.85*a;
%     b = 0.98*b;
% end
end