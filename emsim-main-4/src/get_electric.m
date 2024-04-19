function IV = get_electric( lam, T )
%% FUNCTION ATR1D 
% calculates absorptance, transmittance and reflectance for 1D multilayered stack
% -------------------------------------------------------------------------
%% INPUT
% -------------------------------------------------------------------------
% T
% -------------------------------------------------------------------------
%% OUTPUT
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Constants
% -------------------------------------------------------------------------

if size(lam,1) == 1
    lam = lam';
end

nm = 1e-9;
ec = 1.60217663*1e-19;
h = 6.62607015*1e-34;
c = 3e8;
k = 1.380649*1e-23;
T_room = 300;
n = 1;
I0 = 55*1e-12;
cjsc = ec/(h*c);
cvoc = n*k*T_room/ec;

load("AM15.mat");
AM15_lam = interp1(AM15{"wavelength"},AM15{"direct"},lam);
load("IQE.mat");
IQE_lam = interp1(IQE(:,1),IQE(:,2),lam);

IV.jsc = trapz(lam*nm,cjsc*AM15_lam.*IQE_lam.*T.*lam/10);
IV.Voc = cvoc*log(IV.jsc/I0+1);
IV.V = linspace(0, IV.Voc,1000);
IV.I = IV.jsc - I0*(exp(IV.V/cvoc)-1);
IV.P = IV.I.*IV.V;
IV.Pmax = max(IV.P);
% -------------------------------------------------------------------------
end