%% claculation of CO2 FLUX

% calculate pCO2 using DIC and TA measurements
% import the calcualted DIC and TA
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat')

Kusu.silicate (isnan(Kusu.silicate)) = 5;
Kusu.phosphate (isnan(Kusu.phosphate)) = 0.13;

%% Calculate the pCO2 using DIC and TA
par1type =    1; % The first parameter supplied is of type "1", which is "alkalinity"
par1     =   Kusu.TA; % value of the first parameter
par2type =    2; % The 2nd parameter supplied is of type "2", which is "DIC"
par2     =   Kusu.DIC; % value of the second parameter, which is a long vector of different DIC's!
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
A=CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);

Kusu.pCO2_calc = A(:,4);

%%
% F = k660 K0 (CO2_water - CO2_air)
% pCO2_air from NOAA ESRL Carbon Cycle Coop- erative Global Air Sampling Network (Dlugokencky et al., 2018)
% k660 from Wanninkhof et al. 2009
% K0 from Weiss 1974, lnK0 = A1+A2*(100/T) + A3 * ln(T/100) + S* (B1 + B2/(T/100) + B3*(T/100)^2)
% T is absolute temperature

% function format: K0 = sol_calc (T, Sal, mode)         k660 (cm/hr) = k_calc (T,
% wind speed (unit:m/s));

T = Kusu.temperature;  %temperature, unit:Celcius
salinity = Kusu.salinity;
u = Kusu.wind_speed;  % wind speed, unit: m/s, or we can use the value from each day
pCO2_w = Kusu.pCO2_calc .* 10^-6;
atmCO2 = mean([405.57 407.21 407.10 406.07 407.46 407.90 405.16 403.39 405.37 407.15 407.24 409.57]); %the average of Bukit Kototabang in 2019
pCO2_a = atmCO2 * (10^-6);   %unit:atm, pCO2 in the air, Wit et al. 2018

K0 = sol_calc (T, salinity, 2);  %calculate the solubility constant, input temperature, salinity and unit (1 means mol/kg.atm; 2 means mol/L.atm)
k660 = k_calc (T, u);    %calculate the gas transfer velocity using schmidt number, input temperature and wind speed

% calculate the CO2 flux for each day whenever we have pCO2 data
Kusu.Cflux = k660 .* K0 .* (pCO2_w - pCO2_a) .* 10 .* 24;    %unit: mol/(day*m2)

%% plot
figure('color','w');
plot(Kusu.Date, Kusu.Cflux, 'ko','markerfacecolor','k','markersize',8)
ylabel('CO2 flux (mol/day/m2)')
xlabel('Time')
xlim([datetime('01-Oct-2017') datetime('31-dec-2020')])
ylim([0 0.05])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
set(gca,'xtick',XTICK, 'xticklabel',datestr(XTICK,'mmm-yyyy'),'ytick',[0 0.01 0.02 0.03 0.04 0.05])
set(gca,'fontsize',14)

figure('color','w');
plot(Kusu.Date, Kusu.wind_speed, 'ko','markerfacecolor','k')
ylabel('wind speed(m/s)')
xlabel('Time')
xlim([datetime('01-Oct-2017') datetime('31-dec-2020')])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
set(gca,'xtick',XTICK, 'xticklabel',datestr(XTICK,'mmm-yyyy'))
set(gca,'fontsize',14)


figure('color','w');
plot(Kusu.Date, Kusu.pCO2_calc, 'ko','markerfacecolor','k')
ylabel('pCO2(uatm)')
xlabel('Time')
xlim([datetime('01-Oct-2017') datetime('31-dec-2020')])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
set(gca,'xtick',XTICK, 'xticklabel',datestr(XTICK,'mmm-yyyy'))
set(gca,'fontsize',14)


%% plot the pCO2 data with logger data
% import the logger data
load('F:\NTU\Research\Kusu Hantu Biogeochem\SAMI pCO2 sensor\processed data\pCO2 dailymean\dailymean_co2','dailymean_co2')

%%
figure('color','w');
plot(dailymean_co2.date, dailymean_co2.co2, 'o','color',[102 179 255]/255, 'markerfacecolor', [102 179 255]/255,'markersize',8);hold on
plot(Kusu.Date, Kusu.pCO2_calc, 'o','color',[255 124 124]/255,'markerfacecolor',[255 124 124]/255,'markersize',8);
ylabel('pCO2(uatm)')
xlabel('Time')
xlim([datetime('01-Oct-2017') datetime('31-dec-2020')])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
set(gca,'xtick',XTICK, 'xticklabel',datestr(XTICK,'mmm-yyyy'))
set(gca,'fontsize',14)

%% calculate the annual CO2 flux
%for 2019
startdate = datetime('1-jan-2019');
enddate = datetime('31-dec-2019');

temp1 = Kusu.Cflux (Kusu.Date > startdate & Kusu.Date < enddate & ~isnan(Kusu.Cflux));
temp2 = Kusu.Date(Kusu.Date > startdate & Kusu.Date < enddate & ~isnan(Kusu.Cflux));

joey = table();
joey.DT = vertcat(startdate,temp2,enddate);
joey.Cflux = vertcat(mean(Kusu.Cflux([28 29])),temp1,mean(Kusu.Cflux([53 54]))); % calculate the mean 

joey.datenum = round(datenum(joey.DT - startdate));
figure;plot(joey.DT, joey.Cflux, 'k-')
co2_flux_2019 = trapz(joey.datenum, joey.Cflux);  % mol/m2/yr

%% for 2020
startdate = datetime('1-jan-2020');
enddate = datetime('31-dec-2020');

temp1 = Kusu.Cflux (Kusu.Date > startdate & Kusu.Date < enddate & ~isnan(Kusu.Cflux));
temp2 = Kusu.Date(Kusu.Date > startdate & Kusu.Date < enddate & ~isnan(Kusu.Cflux));

joey = table();
joey.DT = vertcat(startdate,temp2,enddate);
joey.Cflux = vertcat(mean(Kusu.Cflux([53 54])),temp1,mean(Kusu.Cflux([66 67]))); % calculate the mean 

joey.datenum = round(datenum(joey.DT - startdate));
figure;plot(joey.DT, joey.Cflux, 'k-')
co2_flux_2020 = trapz(joey.datenum, joey.Cflux);  % mol/m2/yr

%% calculate the flux during the SW monsoon season
startdate = datetime('28-may-2019');
enddate = datetime('26-sep-2019');

subjoey = fulldata_co2 (fulldata_co2.date > startdate & fulldata_co2.date < enddate & ~isnan(fulldata_co2.flux),:);
subjoey.datenum = round(datenum( subjoey.date - startdate));
figure;plot(subjoey.date, subjoey.flux, 'k-')
co2_flux_sw = trapz(subjoey.datenum, subjoey.flux);  % mol/m2/yr
co2_flux_areatotal_sw = co2_flux_sw * 12 * 127674 * (1000)^2 / (10)^12; % Tg/yr

% annual flux: 4.31 mol/m2/yr,  6.59 TgC/yr
% The CO2 flux during the SW monsoon season accounts for 62% of the annual
% flux

%% monte carlo simulation to estimate the uncertainty
rep = 10000;
MC = struct();
m = size(Kusu,1);
% STEP 1: generate a repition matrix for the parameters
% generate a sereis of random numbers for each day, dont apply the same random
% number series to all the dates

Kusu.k660 = k_calc (Kusu.temperature, Kusu.wind_speed);

for n = 1:m
    MC.salinity(n,1:rep) = randn(1,rep) .* 0.01 + Kusu.salinity(n);
    MC.temperature(n,1:rep) = randn(1,rep) .* 0 + Kusu.temperature(n);
    MC.DIC (n,1:rep) = randn(1,rep) .* Kusu.DIC(n) .* 0.0015 + Kusu.DIC(n);
    MC.TA (n,1:rep) = randn(1,rep) .* Kusu.TA(n) .* 0.0013 + Kusu.TA(n);
    MC.silicate (n,1:rep) = randn(1,rep) .* 0 + Kusu.silicate(n);
    MC.phosphate (n,1:rep) = randn(1,rep) .* 0 + Kusu.phosphate(n);
    MC.k660(n, 1:rep) = randn(1,rep) .* Kusu.k660(n) .* 0.2 + Kusu.k660(n);  % the gas transfer velocity has 20% uncertainty
end


%% STEP 2: repeat the calculation by CO2SYS to calculate the pCO2
% each column corresponds to each trial
for n = 1:rep
    par1type =    1; % The first parameter supplied is of type "1", which is "alkalinity"
    par1     =   MC.TA(:,n); % value of the first parameter
    par2type =    2; % The 2nd parameter supplied is of type "2", which is "DIC"
    par2     =   MC.DIC(:,n); % value of the second parameter, which is a long vector of different DIC's!
    sal      =   MC.salinity(:,n); % Salinity of the sample
    tempin   =   MC.temperature(:,n); % Temperature at input conditions
    presin   =    5; % Pressure    at input conditions
    tempout  =    MC.temperature(:,n); % Temperature at output conditions - doesn't matter in this example
    presout  =    5; % Pressure    at output conditions - doesn't matter in this example
    sil      =    unitqlo(MC.silicate(:,n), MC.salinity(:,n), MC.temperature(:,n), 1); % Concentration of silicate  in the sample (in umol/kg)
    po4      =    unitqlo(MC.phosphate(:,n), MC.salinity(:,n), MC.temperature(:,n), 1); % Concentration of phosphate in the sample (in umol/kg)
    pHscale  =    1; % pH scale at which the input pH is reported ("1" means "Total Scale")  - doesn't matter in this example
    k1k2c    =    4; % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
    kso4c    =    1; % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")
    
    % Do the calculation. See CO2SYS's help for syntax and output format
    A=CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);
    MC.pCO2_calc(:,n) = A(:,4);
end

%% STEP 3: repeat the calculation for the CO2 flux
% %generate a random number series for the wind speed
% for n = 1:m
%     MC.u(n,1:rep) = randn(1,rep) .* 1.4 + Kusu.wind_speed(n); % the uncertainty for the wind speed is 1.4m/s (Ruf et al., 2019)
% end


MC.pCO2_w = MC.pCO2_calc .* 10^-6;
pCO2_a = atmCO2 * (10^-6);   %unit:atm, pCO2 in the air, Wit et al. 2018

MC.K0 = sol_calc (MC.temperature, MC.salinity, 2);  %calculate the solubility constant, input temperature, salinity and unit (1 means mol/kg.atm; 2 means mol/L.atm)
%MC.k660 = k_calc (MC.temperature, MC.u);    %calculate the gas transfer velocity using schmidt number, input temperature and wind speed

% calculate the CO2 flux for each day whenever we have pCO2 data
MC.Cflux = MC.k660 .* MC.K0 .* (MC.pCO2_w - pCO2_a) .* 10 .* 24;    %unit: mol/(day*m2)


%% repeat the integral for the annual CO2 flux and calculate the standard deviation
% for 2019 data
MC.datenum = datenum(vertcat(datetime('1-jan-2019'),Kusu.Date([29:30 33:53]),datetime('31-dec-2019')) - datetime('01-Jan-2019'));
%do not include the NaNs

for n = 1:rep
    flux_startofyear = mean(MC.Cflux([28 29],n));
    flux_endofyear = mean(MC.Cflux([53 54],n));
    data = vertcat(flux_startofyear, MC.Cflux([29:30 33:53],n), flux_endofyear);
    MC.annual_2019(n) = trapz(MC.datenum, data);
end
MC.annual_2019_std = std(MC.annual_2019);


%%
% for 2020 data
MC.datenum = datenum(vertcat(datetime('1-jan-2020'),Kusu.Date(54:66),datetime('31-dec-2020')) - datetime('01-Jan-2020'));
%do not include the NaNs

for n = 1:rep
    flux_startofyear = mean(MC.Cflux([53 54],n));
    flux_endofyear = mean(MC.Cflux([66 67],n));
    data = vertcat(flux_startofyear, MC.Cflux(54:66,n), flux_endofyear);
    MC.annual_2020(n) = trapz(MC.datenum, data);
end
MC.annual_2020_std = std(MC.annual_2020);

%%
MC.Cflux_std = std(MC.Cflux, 0 ,2);

%%
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat','Kusu');
save('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\pCO2\fulldata_co2.mat');


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



function [K0] = sol_calc(T,S,unit)
    if unit == 1   % unit of K0: mol/kg.atm
        A1 = -60.2409;
        A2 = 93.4517;
        A3 = 23.3585;
        B1 = 0.023517;
        B2 = -0.023656;
        B3 = 0.0047036;
        
     
    elseif unit == 2   %unit of K0: mol/(L.atm)
        A1 = -58.0931;
        A2 = 90.5069;
        A3 = 22.2940;
        B1 = 0.027766;
        B2 = -0.025888;
        B3 = 0.0050578;
    end

    T = T+273.15;
    
    lnK0 = A1+A2.*(100./T) + A3 .* log(T./100) + S .* (B1 + B2*(T./100) + B3.*(T./100).^2);
    K0 = exp(lnK0);
    
end

function [k660] = k_calc(T, u)   %T is in degree Celcius   (Wanninkhof, 2014)
    A = 2116.8;
    B = -136.25;
    C = 4.7353;
    D = -0.092307;
    E = 0.0007555;
    Sc = A + B*T + C*T.^2 + D*T.^3 + E * T .^4;  %Schmidt number dependent on temperature
    k660 = 0.251 .* (u.^2) .* ((Sc/660).^(-0.5));
end