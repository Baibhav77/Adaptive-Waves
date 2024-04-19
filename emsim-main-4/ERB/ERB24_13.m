%% SCRIPT FOR TESTING LIGHT PROPAGATION THROUGH BRAGG WITH GLASS
% -------------------------------------------------------------------------
%% INPUT DATA
% -------------------------------------------------------------------------
lamIQE = load('IQE.mat').lam;
IQE = load('IQE.mat').IQE;
lam = linspace(350, 1050, 851);
IQE = interp1(lamIQE,IQE,lam);
%
lamEQE = load('EQE.mat').lam;
EQE = load('EQE.mat').ref/100;
lam = linspace(350, 1050, 851);
EQE = interp1(lamEQE,EQE,lam);
%
Bragg = repmat(["SiO2", "TiO2"],1,4);
%
dBragg =  [110.18, 25.18, 27.09, 50, 22.83, 33.65, 44.78, 12.68]; % SiO2/TiO2 x4
dgls = 3e6; % thickness of glass 50nm
incoh=9;
%
theta = 0; % angle of incidence, normal
% -------------------------------------------------------------------------
%% CALCULATING SPECTRA
% -------------------------------------------------------------------------
% optimized Bragg put just on top of Si cell with ARC optimized to air
% sun -> air|(SiO2|TiO2)x4|glass|(TiO2|SiO2)x4|air|air-ARC|Si
% -------------------------------------------------------------------------
stack = ["air",Bragg,"gls",fliplr(Bragg),"AZO","air","const=2.0661","aSi"];
n = nstack(stack,lam);
dair = 500; dARC = 72.6; dAZO = 5;
dtst = [dBragg, dgls, fliplr(dBragg),dAZO,dair,dARC]';
[A, T, R] = ATR1D( n, lam, dtst, theta, incoh);
% -------------------------------------------------------------------------
%% PLOTTING RESULTS
% -------------------------------------------------------------------------
EQEPSC = T'.*IQE;
figure();
plot(lam,EQE,lam,EQEPSC,'LineWidth',4);
legend('EQE w glass only','EQE w PSC Fraunhofer', ...
    'Location','south','FontSize',16);
set(gca,'FontSize',14)
xlim([lam(1) lam(end)]); ylim([0 1]);
xlabel('wavelength, nm','FontSize',16)
ylabel('EQE','FontSize',16)
% -------------------------------------------------------------------------
% electrical stuff
% -------------------------------------------------------------------------
nm = 1e-9;
ec = 1.60217663*1e-19;
h = 6.62607015*1e-34;
c = 3e8;
k = 1.380649*1e-23;
T = 300;
n = 1;
I0 = 55*1e-12;
cjsc = ec/(h*c);
cvoc = n*k*T/ec;
load("AM15.mat");
spc = interp1(AM15(:,1),AM15(:,2),lam);

jscgls = trapz(lam*nm,cjsc*spc.*EQE.*lam/10);
Vocgls = cvoc*log(jscgls/I0);
V = linspace(0, Vocgls,1000);
I = jscgls - I0*(exp(V/cvoc)-1);
P = I.*V;
Pmaxgls = max(P);

jsc = trapz(lam*nm,cjsc*spc.*EQEPSC.*lam/10);
Voc = cvoc*log(jsc/I0);
V = linspace(0, Voc,1000);
I = jsc - I0*(exp(V/cvoc)-1);
P = I.*V;
Pmax = max(P);

Pgain = (Pmax/Pmaxgls - 1)*100;

title(sprintf('gain in Pmax = %0.1f percent',Pgain));