% calculate the uncertainty of pH and pCO2
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\DOC mineralizatoin\pH_Data.mat','Kusu','Hantu');

%%
[result, headers, units] = test_errors(Kusu);

%%
Kusu.pH_calc_u = result(:,1);
Kusu.pCO2_calc_u = result(:,2);
Kusu.omega_cal_u = result(:,7);
Kusu.omega_arag_u = result(:,8);

%%
%%
[result, headers, units] = test_errors(Hantu);

%%
Hantu.pH_calc_u = result(:,1);
Hantu.pCO2_calc_u = result(:,2);
Hantu.omega_cal_u = result(:,7);
Hantu.omega_arag_u = result(:,8);


%%
save('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\DOC mineralizatoin\KusuHantuData.mat','Kusu','Hantu');
save('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\DOC mineralizatoin\CO2SYS_uncertainty.mat');

%%
function [E,headers,units] = test_errors (Kusu)

PAR1 = Kusu.TA;    % Alk
PAR2 = Kusu.DIC;    % DIC
PAR1TYPE = 1;
PAR2TYPE = 2;
SAL = Kusu.salinity;
TEMPIN = Kusu.temperature;
TEMPOUT = Kusu.temperature;
PRESIN = 5;
PRESOUT = PRESIN;
SI = Kusu.silicate;
PO4 = Kusu.phosphate;
pHSCALEIN = 1;   % total scale
K1K2CONSTANTS = 4; %  Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
KSO4CONSTANTS = 1;  % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")

ePAR1 = Kusu.TA .* 0.0013;    % Alk
ePAR2 = Kusu.DIC .* 0.0015;    % DIC
eSAL = 0.01;
eTEMP = 0;
eSI = 0;
ePO4 = 0;

% No correlation between PAR1 and PAR2
r=0;
% With default errors on Ks
epK = '';
eBt = '';
[E, headers, units] = errors (PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,...
    ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,epK,eBt,r,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);
end
