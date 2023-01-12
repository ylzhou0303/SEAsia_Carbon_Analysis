%% import data
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat','Kusu','Hantu');
load ('F:\NTU\Research\Kusu Hantu Biogeochem\Seabird Sal\processed data for pH correction\daily mean salinity\dailymean_sal.mat','dm_sal');

clear global model
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\endmember.mat')
global model


%% calculate the pH, omega using the measured DIC and TA
Kusu.silicate (isnan(Kusu.silicate)) = 5;
Kusu.phosphate (isnan(Kusu.phosphate)) = 0.13;


par1type =    1; % The first parameter supplied is of type "1", which is "alkalinity"
par1     =    Kusu.TA; % value of the first parameter
par2type =    2; % The 2nd parameter supplied is of type "2", which is "DIC"
par2     =    Kusu.DIC; % value of the second parameter, which is a long vector of different DIC's!
sal      =    Kusu.salinity; % Salinity of the sample
tempin   =    Kusu.temperature; % Temperature at input conditions
presin   =    5; % Pressure    at input conditions
tempout  =    Kusu.temperature; % Temperature at output conditions - doesn't matter in this example
presout  =    5; % Pressure    at output conditions - doesn't matter in this example
sil      =    unitqlo(Kusu.silicate, Kusu.salinity, Kusu.temperature, 1); % Concentration of silicate  in the sample (in umol/kg)
po4      =    unitqlo(Kusu.phosphate, Kusu.salinity, Kusu.temperature, 1); % Concentration of phosphate in the sample (in umol/kg)
pHscale  =    1; % pH scale at which the input pH is reported ("1" means "Total Scale")  - doesn't matter in this example
k1k2c    =    4; % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    =    1; % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")

% Do the calculation. See CO2SYS's help for syntax and output format
A=CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);

Kusu.pH_calc = A(:,3);
Kusu.pCO2_calc = A(:,4);
Kusu.omiga_cal = A(:,15);
Kusu.omiga_arag = A(:,16);

%%
Hantu.silicate (isnan(Hantu.silicate)) = 5;
Hantu.phosphate (isnan(Hantu.phosphate)) = 0.13;
 
 
par1type =    1; % The first parameter supplied is of type "1", which is "alkalinity"
par1     = Hantu.TA; % value of the first parameter
par2type =    2; % The 2nd parameter supplied is of type "2", which is "DIC"
par2     = Hantu.DIC; % value of the second parameter, which is a long vector of different DIC's!
sal      =   Hantu.salinity; % Salinity of the sample
tempin   =   Hantu.temperature; % Temperature at input conditions
presin   =    5; % Pressure    at input conditions
tempout  =    Hantu.temperature; % Temperature at output conditions - doesn't matter in this example
presout  =    5; % Pressure    at output conditions - doesn't matter in this example
sil      =   unitqlo(Hantu.silicate, Hantu.salinity, Hantu.temperature, 1); % Concentration of silicate  in the sample (in umol/kg)
po4      =    unitqlo(Hantu.phosphate, Hantu.salinity, Hantu.temperature, 1); % Concentration of phosphate in the sample (in umol/kg)
pHscale  =    1; % pH scale at which the input pH is reported ("1" means "Total Scale")  - doesn't matter in this example
k1k2c    =    4; % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    =    1; % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")
 
% Do the calculation. See CO2SYS's help for syntax and output format
A=CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);
 
Hantu.pH_calc = A(:,3);
Hantu.pCO2_calc = A(:,4);
Hantu.omiga_cal = A(:,15);
Hantu.omiga_arag = A(:,16);

%% calculate the mean and SD for the Kusu Hantu data of omega during late NE Monsoon and intermonsoon
temp = vertcat(Kusu(:,[1 2 45 49]), Hantu(:, [1 2 36:37]));
temp_mean = mean(temp.omiga_cal([33:36 56:57 97:100 117:118]))
temp_sd = std((temp.omiga_cal([33:36 56:57 97:100 117:118])))
%% CONSERVATIVE MIXING PH
% Here I use the salinity data from the seabird CTD logger to generate a dataset with conservative mixing DIC 
% and TA with a daily resolution, so that I can have a conservative mixing
% pH dataset

dm_sal.salinity(807:837) = NaN;  % these salinity data are believed to be wrong due to large discrepancy with the Valeport sensor (04Nov2018 to 04Dec2018)



%% calculate the conservative mixing DIC and TA, using the seabird CTD data, and then calculate the conservative mixing pH
dm_sal.DIC_cons = DIC_cons_calc(dm_sal.salinity);
dm_sal.TA_cons = TA_cons_calc(dm_sal.salinity);

par1type =    1; % The first parameter supplied is of type "1", which is "alkalinity"
par1     =   dm_sal.TA_cons; % value of the first parameter
par2type =    2; % The 2nd parameter supplied is of type "2", which is "DIC"
par2     =   dm_sal.DIC_cons; % value of the second parameter, which is a long vector of different DIC's!
sal      =   dm_sal.salinity; % Salinity of the sample
tempin   =   dm_sal.temperature; % Temperature at input conditions
presin   =    5; % Pressure    at input conditions
tempout  =    dm_sal.temperature; % Temperature at output conditions - doesn't matter in this example
presout  =    5; % Pressure    at output conditions - doesn't matter in this example
sil      =    5; % Concentration of silicate  in the sample (in umol/kg), mean of the dataset used as a conservative value
po4      =    0.13; % Concentration of phosphate in the sample (in umol/kg), mean of the dataset used as a conservative value
pHscale  =    1; % pH scale at which the input pH is reported ("1" means "Total Scale")  - doesn't matter in this example
k1k2c    =    4; % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    =    1; % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")

% Do the calculation. See CO2SYS's help for syntax and output format
A = CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);
dm_sal.pH_cons = A(:,3);
dm_sal.omiga_cal = A(:,15);
dm_sal.omiga_arag = A(:,16);

%% Also calculate the conservative mixing pH using the Valeport CTD salinity, because for some dates, the seabird logger did not work properly
Kusu.DIC_cons = DIC_cons_calc(Kusu.salinity);
Kusu.TA_cons = TA_cons_calc(Kusu.salinity);

par1type =    1; % The first parameter supplied is of type "1", which is "alkalinity"
par1     =   Kusu.TA_cons; % value of the first parameter
par2type =    2; % The 2nd parameter supplied is of type "2", which is "DIC"
par2     =   Kusu.DIC_cons; % value of the second parameter, which is a long vector of different DIC's!
sal      =   Kusu.salinity; % Salinity of the sample
tempin   =   Kusu.temperature; % Temperature at input conditions
presin   =    5; % Pressure    at input conditions
tempout  =    Kusu.temperature; % Temperature at output conditions - doesn't matter in this example
presout  =    5; % Pressure    at output conditions - doesn't matter in this example
sil      =   unitqlo(Kusu.silicate, Kusu.salinity, Kusu.temperature, 1); % Concentration of silicate  in the sample (in umol/kg)
po4      =    unitqlo(Kusu.phosphate, Kusu.salinity, Kusu.temperature, 1); % Concentration of phosphate in the sample (in umol/kg)
pHscale  =    1; % pH scale at which the input pH is reported ("1" means "Total Scale")  - doesn't matter in this example
k1k2c    =    4; % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    =    1; % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")

% Do the calculation. See CO2SYS's help for syntax and output format
A2 = CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);
Kusu.pH_cons = A2(:,3);
Kusu.omiga_cal_cons = A2(:,15);
Kusu.omiga_arag_cons = A2(:,16);


%% load the measured dailymean pH from the logger data with quality control
load('F:\NTU\Research\Kusu Hantu Biogeochem\SeaFET pH\pH New reprocessed data\dailymean pH\dm_ph_with_QC.mat')

% remove the data that I think is wrong (although I did detailed quality control for the dataset, I still found these dates with suspicious data)
dm_ph{[845:855 221:247 1016:1024 197 1124:1128],2:3} = NaN;
dm_ph{[772:838],2:3} = NaN;

%% PLOT
%conservatie mixing pH, compile the logger data with the Kusu calculated data
ph_cons = table();
ph_cons.DT = vertcat(dm_sal.DT(~isnan(dm_sal.pH_cons)), Kusu.Date(~isnan(Kusu.pH_cons)));
ph_cons.ph = vertcat(dm_sal.pH_cons(~isnan(dm_sal.pH_cons)), Kusu.pH_cons(~isnan(Kusu.pH_cons)));
ph_cons.omiga_cal = vertcat(dm_sal.omiga_cal(~isnan(dm_sal.pH_cons)), Kusu.omiga_cal_cons(~isnan(Kusu.pH_cons)));
ph_cons.omiga_arag = vertcat(dm_sal.omiga_arag(~isnan(dm_sal.pH_cons)), Kusu.omiga_arag_cons(~isnan(Kusu.pH_cons)));
ph_cons = sortrows(ph_cons,'DT');


%% choose which sensor to plot
dm_ph.compiled = vertcat(dm_ph.ext(1:771),dm_ph.int(772:923),dm_ph.ext(924:end));


%%
figure('color','w');
col_cons = [97 227 85]/255;
col_seafet = [102 179 255]/255;
col_calc = [255 124 124]/255;


h3 = plot(ph_cons.DT, ph_cons.ph,'o','markeredgecolor',col_cons,'markerfacecolor',col_cons,'markersize',3); hold on

% h3 = plot(dm_sal.DT, dm_sal.pH_cons, 'go', 'markerfacecolor', 'g');hold on 
% plot(Kusu.Date, Kusu.pH_cons, 'go', 'markerfacecolor', 'g');hold on

% select the more reliable pH sensor for each time point
h1 = plot(dm_ph.DT(1:771),dm_ph.ext(1:771),'o','markerfacecolor',col_seafet,'markeredgecolor',col_seafet ,'markersize',5);   %use external pH until row 771, which is Nov 6th 2018
hold on
plot(dm_ph.DT(839:923),dm_ph.int(839:923),'o','markerfacecolor',col_seafet,'markeredgecolor',col_seafet,'markersize',5);   % use the internal pH since Nov 7th 2018 
hold on
plot(dm_ph.DT(924:end),dm_ph.ext(924:end),'bo','markerfacecolor',col_seafet,'markeredgecolor',col_seafet ,'markersize',5);   % use the external pH since Aug 2nd 2019 
hold on
h2=plot(datetime(datestr(Kusu.Date,0)), Kusu.pH_calc,'o','markeredgecolor',col_calc,'markerfacecolor',col_calc,'markersize',5); hold on   %calculated pH


xlabel('Time');
ylabel('pH')
%title('The pH record from Kusu Island, Singapore (daily mean)')
set(gca,'fontsize',14);

% xlim([datetime('01-oct-2017') datetime('31-Dec-2020')])
% ylim([7.80 8.05])
% XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
% YTICK = [7.8 7.85 7.9  7.95  8  8.05];
% set(gca,'xtick',datetime(XTICK), 'xticklabel',datestr(XTICK,'mmm-yy'),'ytick',YTICK)
% legend ([h1 h2 h3],{'SeaFET','TA+DIC','mixing model'})


% for full dataset
xlim([datetime('01-jul-2015') datetime('31-Dec-2020')])
ylim([7.80 8.05])
XTICK = datetime({'15-jul-2015','15-jan-2016','15-jul-2016','15-jan-2017','15-jul-2017','15-jan-2018','15-jul-2018','15-jan-2019','15-jul-2019','15-jan-2020','15-jul-2020'});
YTICK = [7.8  7.85  7.9  7.95  8  8.05];
set(gca,'xtick',datetime(XTICK), 'xticklabel',datestr(XTICK,'mmm-yy'),'ytick',YTICK)
legend ([h1 h2 h3],{'SeaFET','TA+DIC','mixing model'})

%% plot the omega
% omega of calcite
figure('color','w');
% col_cons = [143	207	209]/255;
% col_meas = [238 187 77]/255;
col_cons = [97 227 85]/255;
col_meas = [102 179 255]/255;


% conservative mixing
plot(ph_cons.DT, ph_cons.omiga_cal, 'o','markeredgecolor',col_cons,'markerfacecolor',col_cons,'markersize',6);hold on
%actual measured data
plot(Kusu.Date, Kusu.omiga_cal,'o','markeredgecolor',col_meas,'markerfacecolor',col_meas,'markersize',8);hold on
plot(Hantu.Date, Hantu.omiga_cal,'o','markeredgecolor',col_meas,'markerfacecolor',col_meas,'markersize',8);

xlim([datetime('01-Oct-2017') datetime('31-Dec-2020')])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
ylim([3 5])
YTICK = [3 3.4 3.8 4.2 4.6 5];
set(gca,'xtick',datetime(XTICK), 'xticklabel',datestr(XTICK,'mmm-yy'),'YTICK',YTICK)
set(gca,'fontsize',14)

%% omega of aragonite
figure('color','w');
% col_cons = [117 207 184]/255;
% col_meas = [246	171	108]/255;
col_cons = [97 227 85]/255;
col_meas = [102 179 255]/255;


plot(ph_cons.DT, ph_cons.omiga_arag, 'o','markeredgecolor',col_cons,'markerfacecolor',col_cons,'markersize',6);hold on
plot(Kusu.Date, Kusu.omiga_arag,'o','markeredgecolor',col_meas,'markerfacecolor',col_meas,'markersize',8);hold on
plot(Hantu.Date, Hantu.omiga_arag,'o','markeredgecolor',col_meas,'markerfacecolor',col_meas,'markersize',8);

xlim([datetime('01-Oct-2017') datetime('31-Dec-2020')])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
ylim([1.8 3.4])
YTICK = [1.8 2.2 2.6 3.0 3.4];
set(gca,'xtick',datetime(XTICK), 'xticklabel',datestr(XTICK,'mmm-yy'),'YTICK',YTICK)
set(gca,'fontsize',14)


%%
save('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\DOC mineralizatoin\pH_Data.mat');
save('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData','Kusu','Hantu');



%% nested functions
function [converted] = unitqlo(original,E2,F2, mode)  %E2 is salinity, F2 is temperature
    %convert the DIC from mol/kg to mol/L, to fit the calculation of CO2
    %flux, which is in the unit of mol/L
    density=(999.842594+0.06793952*F2-0.00909529*F2.^2+0.0001001685*F2.^3-0.000001120083*F2.^4+0.000000006536332*F2.^5+(0.824493-0.0040899*F2+0.000076438*F2.^2-0.00000082467*F2.^3+0.0000000053875*F2.^4).*E2+(-0.00572466+0.00010227*F2-0.0000016546*F2.^2).*E2.^1.5+0.00048314*E2.^2)/1000;
    
    if mode == 1   %Mode 1 is to convert from umol/L to umol/kg
        converted = original ./ density;
    elseif mode == 2 %Mode 2 is to convert from umol/kg to umol/L
        converted = original.*density;
    end
end



function DIC_cons = DIC_cons_calc(salinity)
    global model
    f_river = 1 - salinity ./ model.sal_mar;
    f_marine = salinity ./ model.sal_mar;
    DIC_cons = f_river .* model.dic_riv + f_marine .* model.dic_mar;
end
 
function TA_cons = TA_cons_calc (salinity)
    global model
    f_river = 1 - salinity ./ model.sal_mar;
    f_marine = salinity ./ model.sal_mar;
    TA_cons = f_river .* model.ta_riv + f_marine .* model.ta_mar;
end
