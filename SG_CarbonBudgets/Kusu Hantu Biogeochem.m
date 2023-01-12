%%
% Kusu = readtable('F:\NTU\Research\Kusu Hantu Biogeochem\Kusu Hantu Biogeochem summary\Kusu Hantu 5m Biogeochem data summary.xlsx','sheet','Kusu');
% Hantu = readtable('F:\NTU\Research\Kusu Hantu Biogeochem\Kusu Hantu Biogeochem summary\Kusu Hantu 5m Biogeochem data summary.xlsx','sheet','Hantu');
% 
% Kusu = Kusu(strcmp(Kusu.location,'Kusu'),:);
% Hantu = Hantu(strcmp(Hantu.location,'Hantu'),:);

%% import data
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat','Kusu','Hantu');
%% d13C DOC
figure('color','w')
plot(Kusu.Date, Kusu.d13cDOC, 'o', 'markeredgecolor', [102 179 255]/255,'markerfacecolor', [102 179 255]/255, 'markersize', 8);hold on
plot(Hantu.Date, Hantu.d13cDOC, 'o', 'markeredgecolor', [102 179 255]/255,'markerfacecolor', [102 179 255]/255, 'markersize', 8);

ylim([-26 -19.5])
xlim([datetime('01-Oct-2017') datetime('31-Dec-2020')])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
YTICK = [-26 -25 -24 -23 -22 -21 -20];
set(gca,'xtick',XTICK, 'xticklabel',datestr(XTICK,'mmm-yy'),'ytick',YTICK)
set(gca,'fontsize',14)
%% salinity
%Kusu salinity
a=8;
figure ('color','w')
plot(datetime(datestr(Kusu.Date,0)),Kusu.S_Valeport,'bo','markerfacecolor','b','markersize',a)
hold on
plot(datetime(datestr(Hantu.Date,0)), Hantu.S_Valeport, 'ko', 'markerfacecolor', 'k','markersize',a)
ylabel('DOC (umol/L)')
%xlabel('Date')
xlim([datetime('01-Oct-2017') datetime('31-Dec-2020')])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
set(gca,'xtick',datetime(XTICK), 'xticklabel',datestr(XTICK,'mmm-yy'))
title('Kusu Hantu salinity')
legend('Kusu','Hantu')
set(gca,'fontsize',14)


%% DOC
Hantu.DOC([1 6 7]) = NaN;

a=10;
figure ('color','w')
plot(Kusu.Date,Kusu.DOC,'o', 'markeredgecolor', [102 179 255]/255,'markerfacecolor', [102 179 255]/255,'markersize',a)
hold on
plot(Hantu.Date, Hantu.DOC, 'o', 'markeredgecolor', [102 179 255]/255,'markerfacecolor', [102 179 255]/255,'markersize',a);hold on

xlim([datetime('01-Oct-2017') datetime('31-Dec-2020')])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
set(gca,'xtick',datetime(XTICK), 'xticklabel',datestr(XTICK,'mmm-yy'))
set(gca,'fontsize',14)

%% Kusu a350
figure('color','w')
plot(Kusu.Date,Kusu.a350,'o', 'markeredgecolor', [102 179 255]/255,'markerfacecolor', [102 179 255]/255,'markersize',a)
hold on
plot(Hantu.Date,Hantu.a350,'o', 'markeredgecolor', [102 179 255]/255,'markerfacecolor', [102 179 255]/255,'markersize',a)


xlim([datetime('01-Oct-2017') datetime('31-Dec-2020')])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
set(gca,'xtick',datetime(XTICK), 'xticklabel',datestr(XTICK,'mmm-yy'))
set(gca,'fontsize',14)


%% Kusu Hantu s275-295
figure('color','w')
plot(Kusu.Date,Kusu.s275_295,'o', 'markeredgecolor', [102 179 255]/255,'markerfacecolor', [102 179 255]/255,'markersize',a)
hold on
plot(Hantu.Date,Hantu.s275_295,'o', 'markeredgecolor', [102 179 255]/255,'markerfacecolor', [102 179 255]/255,'markersize',a)

YTICK = [0.012 0.016 0.02 0.024 0.028 0.032];
ylim([0.012 0.033])
xlim([datetime('01-Oct-2017') datetime('31-Dec-2020')])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
set(gca,'xtick',datetime(XTICK), 'xticklabel',datestr(XTICK,'mmm-yy'),'ytick',YTICK);
set(gca,'fontsize',14)

%% Kusu Hantu SUVA254
figure('color','w')
plot(Kusu.Date,Kusu.SUVA254,'o', 'markeredgecolor', [102 179 255]/255,'markerfacecolor', [102 179 255]/255,'markersize',a)
hold on
plot(Hantu.Date,Hantu.SUVA254,'o', 'markeredgecolor', [102 179 255]/255,'markerfacecolor', [102 179 255]/255,'markersize',a)

ylabel('SUVA_2_5_4 (L mg^-^1 m^-^1)')
xlim([datetime('01-Oct-2017') datetime('31-Dec-2020')])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
set(gca,'xtick',datetime(XTICK), 'xticklabel',datestr(XTICK,'mmm-yy'))
%title('Kusu Hantu SUVA254')
set(gca,'fontsize',14)


%% Kusu Hantu chl-a
figure('color','w')
plot(Kusu.Date(26:end),Kusu.chl_a(26:end),'o','markeredgecolor',[102 179 255]/255,'markerfacecolor',[102 179 255]/255,'markersize',8)
hold on
plot(Hantu.Date(25:end),Hantu.chl_a(25:end),'o','markeredgecolor',[102 179 255]/255,'markerfacecolor',[102 179 255]/255,'markersize',8)

ylabel('chl-a (ug/L)')
%ylim([0 3.5])
xlim([datetime('01-Oct-2017') datetime('31-Dec-2020')])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
set(gca,'xtick',datetime(XTICK), 'xticklabel',datestr(XTICK,'mmm-yy'))
set(gca,'fontsize',14)

%% Kusu Hantu DIC
figure('color','w')

plot(datetime(datestr(Kusu.Date,0)),Kusu.DIC,'bo','markerfacecolor','b','markersize',a)
hold on
plot(datetime(datestr(Hantu.Date,0)),Hantu.DIC,'ko','markerfacecolor','k','markersize',a)
hold on
%plot(Kusu.Date, Kusu.DIC_calc,'rx','markersize',a)

ylabel('salinity-normalized DIC (umol/kg)')
xlim([datetime('01-Jul-2018') datetime('30-oct-2019')])

Xtick = {'15-Oct-2017','15-Jan-2018','15-Apr-2018','15-Jul-2018','15-Oct-2018','15-Jan-2019','15-Apr-2019','15-Jul-2019','15-oct-2019'}
set(gca,'xtick',datetime(Xtick), 'xticklabel',datestr(Xtick,'mmm-yy'),'fontsize',14)
legend('Kusu','Hantu')
%title('Kusu Hantu DIC')

%% Kusu Hantu TA
figure('color','w')
plot(datetime(datestr(Kusu.Date,0)),Kusu.TA,'bo','markerfacecolor','b','markersize',a)
hold on
plot(datetime(datestr(Hantu.Date,0)),Hantu.TA,'ko','markerfacecolor','k','markersize',a)
hold on
%plot(Kusu.Date,Kusu.TA_calc,'rx','markersize',a)

ylabel('total alkalinity (umol/kg)')
xlim([datetime('01-Jul-2018') datetime('30-oct-2019')])
Xtick = {'15-Oct-2017','15-Jan-2018','15-Apr-2018','15-Jul-2018','15-Oct-2018','15-Jan-2019','15-Apr-2019','15-jul-2019','15-oct-2019'};
set(gca,'xtick',datetime(Xtick), 'xticklabel',datestr(Xtick,'mmm-yy'),'fontsize',14)
legend('Kusu','Hantu')
title('Kusu Hantu total alkalinity')


%% TA vs salinity
figure('color','w')
a=10
h1 = plot(Kusu.salinity,Kusu.TA,'ko','markerfacecolor','k','markersize',a);
hold on
plot(Hantu.salinity,Hantu.TA,'ko','markerfacecolor','k','markersize',a);
hold on
h2 = plot(Kusu.salinity(strcmp(Kusu.season,'NE monsoon')==1),Kusu.TA(strcmp(Kusu.season,'NE monsoon')==1),'bo','markerfacecolor','b','markersize',a);
hold on
plot(Hantu.salinity(strcmp(Hantu.season,'NE monsoon')==1),Hantu.TA(strcmp(Hantu.season,'NE monsoon')==1),'bo','markerfacecolor','b','markersize',a);
hold on
h3 = plot(Kusu.salinity(strcmp(Kusu.season,'SW monsoon')==1),Kusu.TA(strcmp(Kusu.season,'SW monsoon')==1),'go','markerfacecolor','g','markersize',a);
hold on
plot(Hantu.salinity(strcmp(Hantu.season,'SW monsoon')==1),Hantu.TA(strcmp(Hantu.season,'SW monsoon')==1),'go','markerfacecolor','g','markersize',a);
hold on
%plot(Kusu.S_Valeport,Kusu.TA_calc,'rx','markersize',a);


ylim([2000 2200])
legend ([h1 h2 h3],{'inter-monsoon','NE monsoon','SW monsoon'})
set(gca,'fontsize',14)
%% linear regression
% TA = 46.153 * Sal + 631.72
sal_mod = [29.6:0.1:32.7];
plot (sal_mod, 46.153*sal_mod+631.72, 'k-', 'linewidth', 3)

%% Hantu a350
subplot(2,2,2)
plot(Hantu.Date,Hantu.a350,'ko','markerfacecolor','k')
ylabel('a_3_5_0 (m^-^1)')
%xlabel('Date')
xlim([datetime('01-Oct-2017') datetime('15-Dec-2018')])
set(gca,'xtick',datetime({'01-Oct-2017','01-Mar-2018','01-Jul-2018','01-Dec-2018'}))
%set(gca,'xticklabel',datestr({'01-Oct-2017','01-Mar-2018','01-Jul-2018','01-Dec-2018'},'mmm-yyyy'))
title('Hantu a350')
set(gca,'fontsize',14)


%Hantu s275-295
subplot(2,2,4)
plot(Hantu.Date,Hantu.s275_295,'ko','markerfacecolor','k')
ylabel('S_2_7_5_-_2_9_5')
xlim([datetime('01-Oct-2017') datetime('15-Nov-2018')])
xlim([datetime('01-Oct-2017') datetime('15-Dec-2018')])
set(gca,'xtick',datetime({'01-Oct-2017','01-Mar-2018','01-Jul-2018','01-Dec-2018'}))
set(gca,'xticklabel',datestr({'01-Oct-2017','01-Mar-2018','01-Jul-2018','01-Dec-2018'},'mmm-yyyy'))
title('Hantu S275-295')
set(gca,'fontsize',14)


%% loggers
%% seabird dailymean salinity
% import logger data
cd('F:\NTU\Research\Kusu Hantu Biogeochem\Seabird Sal\processed data for pH correction\daily mean salinity')
S=readtable('daily mean salinity.xlsx');
%remove the data since Nov 01 2018
S=vertcat(S(1:803,:),S(838:end,:));


%Arrange the Xtick
Xtick = {'15-Jul-2015','15-Oct-2015','15-Jan-2016','15-Apr-2016','15-Jul-2016','15-Oct-2016','15-Jan-2017','15-Apr-2017','15-Jul-2017','15-Oct-2017','15-Jan-2018','15-Apr-2018','15-Jul-2018','15-Oct-2018','15-Jan-2019','15-Apr-2019','15-Jul-2019'};
Xtick = datetime(Xtick);

figure('color','w')
plot(S.date,S.salinity,'bo','markerfacecolor','b','markersize',a);
hold on
plot(Kusu.Date,Kusu.S_Valeport,'go','markerfacecolor','g','markersize',a);

title('Salinity record from Kusu Island, Singapore')
xlabel('Time')
ylabel('Salinity')
xlim ([datetime('01-Jan-2017') datetime('15-Aug-2019')])
legend('Seabird CTD','Valeport CTD')
set(gca,'xtick',datetime(Xtick),'xticklabel',datestr(Xtick,'mmm-yy'),'fontsize',14);

%% calculated pH using CO2SYS
addpath('F:\Matlab Files\CO2SYS_MATLAB\CO2SYS-MATLAB-master\src');
Kusu.P_kg = unitqlo(Kusu.phosphate, Kusu.salinity, Kusu.temperature, 1);
Kusu.Si_kg = unitqlo(Kusu.silicate, Kusu.salinity, Kusu.temperature, 1);

% fill in the missing (not measured) nutrient data
Kusu.P_kg (isnan(Kusu.P_kg)) = 0.15;
Kusu.Si_kg (isnan(Kusu.Si_kg)) = 5;

par1type =    1; % The first parameter supplied is of type "1", which is "alkalinity"
par1     =   Kusu.TA; % value of the first parameter
par2type =    2; % The 2nd parameter supplied is of type "2", which is "DIC"
par2     =   Kusu.DIC; % value of the second parameter, which is a long vector of different DIC's!
sal      =   Kusu.salinity; % Salinity of the sample
tempin   =   Kusu.temperature; % Temperature at input conditions
presin   =    5; % Pressure    at input conditions
tempout  =    Kusu.temperature; % Temperature at output conditions - doesn't matter in this example
presout  =    5; % Pressure    at output conditions - doesn't matter in this example
sil      =   Kusu.Si_kg; % Concentration of silicate  in the sample (in umol/kg)
po4      =   Kusu.P_kg; % Concentration of phosphate in the sample (in umol/kg)
pHscale  =    1; % pH scale at which the input pH is reported ("1" means "Total Scale")  - doesn't matter in this example
k1k2c    =    4; % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    =    1; % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")

% Do the calculation. See CO2SYS's help for syntax and output format
A=CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);
Kusu.pH_calc = A(:,3);
Kusu.pCO2_calc = A(:,4);

%% Hantu calculated pH
Hantu.P_kg = unitqlo(Hantu.phosphate, Hantu.salinity, Hantu.temperature, 1);
Hantu.Si_kg = unitqlo(Hantu.silicate, Hantu.salinity, Hantu.temperature, 1);
Hantu.P_kg (isnan(Hantu.P_kg)) = 0.15;
Hantu.Si_kg (isnan(Hantu.Si_kg)) = 5;

par1type =    1; % The first parameter supplied is of type "1", which is "alkalinity"
par1     =   Hantu.TA; % value of the first parameter
par2type =    2; % The 2nd parameter supplied is of type "2", which is "DIC"
par2     =   Hantu.DIC; % value of the second parameter, which is a long vector of different DIC's!
sal      =   Hantu.salinity; % Salinity of the sample
tempin   =   Hantu.temperature; % Temperature at input conditions
presin   =    5; % Pressure    at input conditions
tempout  =    Hantu.temperature; % Temperature at output conditions - doesn't matter in this example
presout  =    5; % Pressure    at output conditions - doesn't matter in this example
sil      =   Hantu.Si_kg; % Concentration of silicate  in the sample (in umol/kg)
po4      =   Hantu.P_kg; % Concentration of phosphate in the sample (in umol/kg)
pHscale  =    1; % pH scale at which the input pH is reported ("1" means "Total Scale")  - doesn't matter in this example
k1k2c    =    4; % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    =    1; % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")
 
% Do the calculation. See CO2SYS's help for syntax and output format
A=CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);
Hantu.pH_calc = A(:,3);
Hantu.pCO2_calc = A(:,4);

%% save data
save('F:\NTU\Research\Kusu Hantu Biogeochem\Kusu Hantu Biogeochem summary\KusuHantuBiogeochem','Kusu','Hantu')
%%
h_unitqlo = @unitqlo
%% pH
% import the logger data
cd('F:\NTU\Research\Kusu Hantu Biogeochem\SeaFET pH\pH New reprocessed data\dailymean pH')
pH=readtable('dailymean_pH.xlsx');


% plot daily mean pH
 %arrange the xtick
Xtick = {'15-Jul-2015','15-Oct-2015','15-Jan-2016','15-Apr-2016','15-Jul-2016','15-Oct-2016','15-Jan-2017','15-Apr-2017','15-Jul-2017','15-Oct-2017','15-Jan-2018','15-Apr-2018','15-Jul-2018','15-Oct-2018','15-Jan-2019','15-Apr-2019','15-Jul-2019'};
Xtick = datetime(Xtick);

% we can choose one of the sensor that we believe gives correct number to plot
figure('color','w');
h1=plot(pH.date(1:771),pH.external(1:771),'bo','markerfacecolor','b','markersize',a);   %use external pH until row 771, which is Nov 6th 2018
hold on

% the data between Nov 2018 and May 2019 are wrong
plot(pH.date(860:end),pH.internal(860:end),'bo','markerfacecolor','b','markersize',a);   % use the internal pH since Nov 7th 2018
hold on
h2=plot(Kusu.Date, Kusu.pH_calc,'go','markerfacecolor','g','markersize',a)


xlabel('Time');
ylabel('pH')
title('pH record at Kusu Island, Singapore')
set(gca,'fontsize',16);
xlim([datetime('01-Jan-2017') datetime('15-aug-2019')])
ylim([7.8 8.05])
set(gca, 'XTick', Xtick)
set(gca,'Xticklabel',datestr(Xtick,'mmm-yy'))
legend ([h1 h2],{'SeaFET logger','calculated pH'})


%% pCO2
%import logger data
cd('F:\NTU\Research\Kusu Hantu Biogeochem\SAMI pCO2 sensor\processed data\pCO2 Matlab workspace time series')
load dailymeanpCO2

pCO2 = pCO2_dailymean;

% import the calculated pCO2
cd('F:\NTU\Research\Kusu Hantu Biogeochem\Kusu Hantu Biogeochem summary')
load KusuHantu_Biogeochem

%remove abnormal wrong values
pCO2 = pCO2(pCO2.dailymean<600,:);

Xtick = {'01-Jul-2015','01-Oct-2015','01-Jan-2016','01-Apr-2016','01-Jul-2016','01-Oct-2016','01-Jan-2017','01-Apr-2017','01-Jul-2017','01-Oct-2017','01-Jan-2018','01-Apr-2018','01-Jul-2018','01-Oct-2018','01-Feb-2019'};
Xtick = datetime(Xtick);

figure('color','w')
plot(datetime(pCO2.date),pCO2.dailymean,'bo','markerfacecolor','b','markersize',10);
hold on
plot(Kusu.Date,Kusu.pCO2_calc,'ro','markerfacecolor','r','markersize',10);

title('pCO2 record from Kusu Island, Singapore')
xlabel('Time')
ylabel('pCO2 (uatm)')
xlim ([datetime('01-Oct-2017') datetime('15-Mar-2019')])
legend ('SAMI pCO2 logger','calculated pCO2')
set(gca,'xtick',datetime(Xtick),'xticklabel',datestr(Xtick,'mmm-yyyy'),'fontsize',14);


%% nutrients
cd('F:\NTU\Research\Kusu Hantu Biogeochem\nutrients')
nutrients = readtable ('Kusu Hantu summary of nutrients.xlsx','sheet','Kusu Hantu');

Kusu_n = nutrients(strcmp(nutrients.location,'Kusu')==1,:);
Hantu_n = nutrients(strcmp(nutrients.location,'Hantu')==1,:);

figure('color','w')
plot(Kusu_n.salinity,Kusu_n.nitrate,'bo','markerfacecolor','b','markersize',a)
hold on
plot(Hantu_n.salinity,Hantu_n.nitrate,'ko','markerfacecolor','k','markersize',a)

Xtick = {'15-jul-2014','15-jul-2015','15-jul-2016','15-jul-2017','15-jul-2018','15-jul-2019'};
set(gca,'xtick',datetime(Xtick),'xticklabel',datestr(Xtick,'mmm-yy'))
ylabel('NO3 (umol/L)')
set(gca,'fontsize',14)

%%
figure('color','w')
plot(Kusu_n.salinity,Kusu_n.nitrite,'bo','markerfacecolor','b','markersize',a)
hold on
plot(Hantu_n.salinity,Hantu_n.nitrite,'ko','markerfacecolor','k','markersize',a)

Xtick = {'15-jul-2014','15-jul-2015','15-jul-2016','15-jul-2017','15-jul-2018','15-jul-2019'};
set(gca,'xtick',datetime(Xtick),'xticklabel',datestr(Xtick,'mmm-yy'))
ylabel('NO2 (umol/L)')
set(gca,'fontsize',14)


%%
figure('color','w')
plot(Kusu_n.date,Kusu_n.Ammonia,'bo','markerfacecolor','b','markersize',a)
hold on
plot(Hantu_n.date,Hantu_n.Ammonia,'ko','markerfacecolor','k','markersize',a)

ylabel('NH3 (umol/L)')
ylim([0 1.2])
set(gca,'fontsize',14)


%%
figure('color','w')
plot(Kusu_n.salinity,Kusu_n.phosphate,'bo','markerfacecolor','b','markersize',a)
hold on
plot(Hantu_n.salinity,Hantu_n.phosphate,'ko','markerfacecolor','k','markersize',a)

ylabel('PO4 (umol/L)')
%ylim([0 1.2])
set(gca,'fontsize',14)

%%
figure('color','w')
plot(Kusu_n.date,Kusu_n.N_PRatio,'bo','markerfacecolor','b','markersize',a)
hold on
plot(Hantu_n.date,Hantu_n.N_PRatio,'ko','markerfacecolor','k','markersize',a)
hold on
plot([datetime('15-jul-2014'), datetime('15-jul-2019')], [16 ,16],'k-')

Xtick = {'15-jul-2014','15-jul-2015','15-jul-2016','15-jul-2017','15-jul-2018','15-jul-2019'};
set(gca,'xtick',datetime(Xtick),'xticklabel',datestr(Xtick,'mmm-yy'))

ylabel('DIN:DIP RATIO')
ylim([0 30])
set(gca,'fontsize',14)


%%
figure('color','w')
plot(Kusu_n.salinity,Kusu_n.silicate,'bo','markerfacecolor','b','markersize',a)
hold on
plot(Hantu_n.salinity,Hantu_n.silicate,'ko','markerfacecolor','k','markersize',a)
xlabel('salinity')
ylabel('silicate (umol/L)')
set(gca,'fontsize',14)

%% Kusu Hantu d13c-dic
cd('F:\NTU\Research\Kusu Hantu Biogeochem\13C DIC')
DIC13C = readtable('20190910_d13C_Mastersheet.xlsx');
DIC13C.Date = datetime(DIC13C.SampleCollectionDate);
figure('color','w')
plot(DIC13C.Date, DIC13C.d13C, 'ko', 'markerfacecolor','k');

ylabel('d13C-DIC')
set(gca,'fontsize',14)
%% Kusu Hantu DIC vs Talk
figure('color','w')

%southwest monsoon season values
DIC_sw = vertcat(Kusu.DIC(strcmp(Kusu.season, 'SW monsoon')==1),Hantu.DIC(strcmp(Hantu.season, 'SW monsoon') == 1));
TA_sw = vertcat(Kusu.TA(strcmp(Kusu.season, 'SW monsoon')==1),Hantu.TA(strcmp(Hantu.season, 'SW monsoon') == 1));
S_sw = vertcat(Kusu.S_Valeport(strcmp(Kusu.season, 'SW monsoon')==1),Hantu.S_Valeport(strcmp(Hantu.season, 'SW monsoon') == 1));

%northeast monsoon season values
DIC_ne = vertcat(Kusu.DIC(strcmp(Kusu.season, 'NE monsoon')==1),Hantu.DIC(strcmp(Hantu.season, 'NE monsoon') == 1));
TA_ne = vertcat(Kusu.TA(strcmp(Kusu.season, 'NE monsoon')==1),Hantu.TA(strcmp(Hantu.season, 'NE monsoon') == 1));
S_ne = vertcat(Kusu.S_Valeport(strcmp(Kusu.season, 'NE monsoon')==1),Hantu.S_Valeport(strcmp(Hantu.season, 'NE monsoon') == 1));


TA_all = vertcat(Kusu.TA, Hantu.TA);
DIC_all = vertcat(Kusu.DIC, Hantu.DIC);
S_all = vertcat(Kusu.S_Valeport, Hantu.S_Valeport);


% normalize the TA to the salinity of 35
TA_sw_n = TA_sw./S_sw*35;
TA_ne_n = TA_ne./S_ne*35;
TA_all_n = TA_all./S_all*35;

%normalize the DIC to the salinity of 35
DIC_sw_n = DIC_sw./S_sw*35;
DIC_ne_n = DIC_ne./S_ne*35;
DIC_all_n = DIC_all./S_all*35;


% % use the normalization approach by Friis 2003
% TA_sw_n = (TA_sw - 298.87)./S_sw*35 + 298.87;
% TA_ne_n = (TA_ne - 298.87)./S_ne*35 + 298.87;
% TA_all_n = (TA_all - 298.87)./S_all*35 + 298.87;
% 
% DIC_sw_n = (DIC_sw - 676.46)./S_sw*35 + 676.46;
% DIC_ne_n = (DIC_ne-676.46)./S_ne*35 + 676.46;
% DIC_all_n = (DIC_all-676.46)./S_all*35 + 676.46;


% now plot the data
figure
plot(DIC_all_n, TA_all_n, 'ko','markerfacecolor','k');
hold on
plot(DIC_sw_n, TA_sw_n, 'go','markerfacecolor','g');
hold on
plot(DIC_ne_n, TA_ne_n, 'bo','markerfacecolor','b');
hold on



% add some reference lines
o_point_x = mean(DIC_all_n,'omitnan');
o_point_y = mean(TA_all_n,'omitnan');
plot(o_point_x, o_point_y, 'ko', 'markersize', 10)
hold on

% CaCO3 formation/dissolution line, delta TA: delta DIC = 2:1
plot([o_point_x o_point_x+30],[o_point_y o_point_y+2*30],'k-')
hold on
plot([o_point_x-30 o_point_x],[o_point_y-2*30 o_point_y],'k-')
hold on

% respiration/photosynthesis: only consider the C:N = 106:16
% respiration of OM will release NO3- so the TA decreases while DIC
% increases
plot([o_point_x o_point_x+30], [o_point_y   o_point_y - 16/106*30], 'k-')   %respiration
hold on
plot([o_point_x-30 o_point_x],[o_point_y + 16/106*30  o_point_y],'k-')   %photosynthesis
hold on

%CO2 gas exchange
plot([o_point_x-30 o_point_x+30],[o_point_y o_point_y],'k-')

xlim([2020 2140])
ylim([2280 2360])
xlabel('nDIC (umol/kg)')
ylabel('nTA (umol/kg)')
%% under different ratio of OM:CaCO3
figure ('color','w')
plot(DIC_all_n, TA_all_n, 'ko','markerfacecolor','k');
hold on
plot(DIC_sw_n, TA_sw_n, 'go','markerfacecolor','g');
hold on
plot(DIC_ne_n, TA_ne_n, 'bo','markerfacecolor','b');
hold on

% add lines
x = [2030 2140]
%plot([2030 2140], [2290 2290+(2140-2030)]  , 'k-');  %delta DIC:delta TA = 1?1
hold on
plot([2030 2140], [2290 2290+0.5*(2140-2030)], 'k-');  %delta DIC:deltaTA = 2:1
hold on
%plot([2030 2140], [2290 2290+2*(2140-2030)], 'k-');  %delta DIC:deltaTA = 1:2
hold on
%plot([2060 2140],[2290 2290+3*(2140-2060)],'k-')  %delta DIC:delta TA = 1:3

ylim([2285 2360])
%% DIC vs salinity
figure('color','w')

%southwest monsoon season values
DIC_sw = vertcat(Kusu.DIC(strcmp(Kusu.season, 'SW monsoon')==1),Hantu.DIC(strcmp(Hantu.season, 'SW monsoon') == 1));
S_sw = vertcat(Kusu.S_Valeport(strcmp(Kusu.season, 'SW monsoon')==1),Hantu.S_Valeport(strcmp(Hantu.season, 'SW monsoon') == 1));

%northeast monsoon season values
DIC_ne = vertcat(Kusu.DIC(strcmp(Kusu.season, 'NE monsoon')==1),Hantu.DIC(strcmp(Hantu.season, 'NE monsoon') == 1));
S_ne = vertcat(Kusu.S_Valeport(strcmp(Kusu.season, 'NE monsoon')==1),Hantu.S_Valeport(strcmp(Hantu.season, 'NE monsoon') == 1));



h1 = plot(vertcat(Kusu.S_Valeport, Hantu.S_Valeport), vertcat(Kusu.DIC, Hantu.DIC), 'ko', 'markerfacecolor','k')
hold on
h2 = plot(S_sw, DIC_sw, 'go', 'markerfacecolor','g')
hold on
h3 = plot(S_ne, DIC_ne, 'bo', 'markerfacecolor', 'b')

legend([h1, h2, h3], {'inter-monsoon', 'SW monsoon', 'NE monsoon'})



%% everything plotted against salinity
figure('color','w')
subplot(3,3,1)
plot(Kusu_raw.S_Valeport, Kusu_raw.DOC, 'bo','markerfacecolor','b')
hold on
plot(Hantu_raw.S_Valeport, Hantu_raw.DOC, 'ko', 'markerfacecolor', 'k')
ylabel('DOC (umol/L)')


subplot(3,3,2)
plot(Kusu_raw.S_Valeport , Kusu_raw.pH_internal__dailyMean , 'bo','markerfacecolor','b')
hold on
ylabel('pH dailymean (internal sensor)')
ylim([7.8 7.95])

subplot(3,3,3)
plot(Kusu_raw.S_Valeport , Kusu_raw.pH_external__dailymean , 'bo','markerfacecolor','b')
hold on
ylabel('pH dailymean (external sensor)')

subplot(3,3,4)
plot(Kusu.salinity, Kusu.pH, 'bo','markerfacecolor','b')
ylabel('pH at Kusu')

subplot(3,3,5)
plot(Kusu_raw.salinity, Kusu_raw.a350, 'bo','markerfacecolor','b')
hold on
plot(Hantu_raw.S_Valeport, Hantu_raw.a350, 'ko', 'markerfacecolor', 'k')
ylabel('a350 (m^-^1)')
ylim([0 2.1])


subplot(3,3,6)
plot(Kusu_raw.salinity, Kusu_raw.SUVA254, 'bo','markerfacecolor','b')
hold on
plot(Hantu_raw.S_Valeport, Hantu_raw.SUVA254, 'ko', 'markerfacecolor', 'k')
ylabel('SUVA254')

subplot(3,3,7)
plot(Kusu_raw.salinity, Kusu_raw.s275_295, 'bo','markerfacecolor','b')
hold on
plot(Hantu_raw.S_Valeport, Hantu_raw.s275_295, 'ko', 'markerfacecolor', 'k')
ylabel('S 275-295')

subplot(3,3,8)
plot(Kusu_raw.salinity, Kusu_raw.s350_400, 'bo','markerfacecolor','b')
hold on
plot(Hantu_raw.S_Valeport, Hantu_raw.s350_400, 'ko', 'markerfacecolor', 'k')
ylabel('S 350-400')

subplot(3,3,9)
plot(Kusu_raw.salinity, Kusu_raw.sloperatio, 'bo','markerfacecolor','b')
hold on
plot(Hantu_raw.S_Valeport, Hantu_raw.sloperatio, 'ko', 'markerfacecolor', 'k')
ylabel('spectral slope ratio')



%% random stuffs
figure;
plot(Kusu.Date, Kusu.pH_calc,'go','markerfacecolor','g');
hold on
plot(Kusu.Date, Kusu.pH_internal,'go','markerfacecolor','b')
xlim([datetime('01-May-2019') datetime('15-Aug-2019')])



%%
figure
m = size(Kusu)
for i = 1:m(1)
    xlim([0.016 0.03])
    ylim([7.8 8.05])
    plot(Kusu.s275_295(1:i), Kusu.pH_compiled(1:i), 'ko','markerfacecolor','k')
    hold on
    plot(Kusu.s275_295(i), Kusu.pH_compiled(i), 'ro','markerfacecolor','r')
    hold on
    fprintf(datestr(Kusu.Date(i)))
    fprintf('\n')
    pause(3)
end



%%
plot(d13c.SampleCollectionDate, d13c.d13C, 'ko', 'markerfacecolor', 'k');hold on
plot(Kusu.Date, Kusu.d13c_pred, 'go', 'markerfacecolor', 'g');hold on
plot(Kusu.Date, Kusu.d13c_mixed, 'bo','markerfacecolor', 'b');hold on

legend ('measured', 'conservative mixing', 'mixed with DOC')

%xlim ([datetime('') )

%%
figure;
plot(Kusu.S_Valeport, Kusu.mDOC, 'ko');hold on
plot(Hantu.S_Valeport, Hantu.mDOC, 'ko')

figure;
plot(Kusu.S_Valeport, Kusu.a350, 'ko');hold on
plot(Hantu.S_Valeport, Hantu.a350, 'ko')

figure;
plot(Kusu.S_Valeport, Kusu.SUVA254, 'ko');hold on
plot(Hantu.S_Valeport, Hantu.SUVA254, 'ko')

figure;
plot(Kusu.Date, Kusu.tDOC./Kusu.DOC, 'ko');hold on
plot(Hantu.Date, Hantu.tDOC./Hantu.DOC, 'ko')

figure;
plot(Kusu.a350, Kusu.tDOC, 'ko');hold on
plot(Hantu.a350, Hantu.tDOC, 'ko')

figure;
plot(Kusu.Date, Kusu.SUVA254, 'ko');hold on
plot(Hantu.Date, Hantu.SUVA254, 'ko')

figure;
plot(Kusu2.Date, Kusu2.d13cDOC, 'ko');hold on
plot(Hantu2.Date, Hantu2.d13cDOC, 'ko')


%%
X = vertcat(Kusu.a350, Hantu.a350);
X = horzcat(ones(size(X)), X);
Y = vertcat(Kusu.tDOC, Hantu.tDOC);

[B,BINT,R,RINT,STATS] = regress(Y,X)
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