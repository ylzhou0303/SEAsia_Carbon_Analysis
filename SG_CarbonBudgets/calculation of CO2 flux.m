%% claculation of CO2 FLUX

% calculate pCO2 using DIC and TA measurements
% import the calcualted DIC and TA
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat')

Kusu.silicate (isnan(Kusu.silicate)) = 5;
Kusu.phosphate (isnan(Kusu.phosphate)) = 0.13;

%% Calculate the pCO2 using DIC and TA
addpath('F:\NTU\Research\Kusu Hantu Biogeochem\CO2SYS_MATLAB\CO2SYS-MATLAB-master\src');
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
clearvars -except Kusu
%% Compile the pCO2 from SAMI logger with the calculated pCO2
% import dailymean pCO2 value measured by the SAMI logger
load('F:\NTU\Research\Kusu Hantu Biogeochem\SAMI pCO2 sensor\processed data\pCO2 dailymean\dailymean_co2','dailymean_co2')

% feed the salinity to the dailymean pCO2 data
load ('F:\NTU\Research\Kusu Hantu Biogeochem\Seabird Sal\processed data for pH correction\daily mean salinity\dailymean_sal.mat','dm_sal');
dm_sal.date = datetime(datestr(dm_sal.date, 'dd-mmm-yyyy'));

dm_sal_hana = dm_sal (dm_sal.date > datetime('19-mar-2018') & dm_sal.date < datetime('14-jun-2018'),:);
m = size(dailymean_co2);
dailymean_co2.sal = zeros(m(1),1);
dailymean_co2.sal(:) = NaN;
for i = 1:m(1)
    temp = dm_sal_hana.salinity(strcmp(cellstr(dm_sal_hana.date),datestr(dailymean_co2.date(i),'dd-mmm-yyyy')) == 1);
    if ~isempty(temp)
        dailymean_co2.sal(i) = temp;
    end
end


%% feed the salinity from the valeport

for i = 1:m(1)
    if isnan(dailymean_co2.sal(i))
        kuk = Kusu.salinity_compiled (strcmp(cellstr(datestr(Kusu.Date, 'dd-mmm-yyyy')), datestr(dailymean_co2.date(i), 'dd-mmm-yyyy') ) == 1);
        if ~isempty(kuk)
            dailymean_co2.sal(i) = kuk;
        end
    end
end


dailymean_co2.sal(isnan(dailymean_co2.sal)) = mean([Kusu.salinity_compiled(8), Kusu.salinity_compiled(10)]);


%% now feed the wind speed into the dailymean pCO2 dataset
m = size(dailymean_co2);

for i = 1:m(1)
    date_dif = abs(datenum(dailymean_co2.date(i) - Kusu.Date));
    dailymean_co2.wind_speed (i) = Kusu.wind_speed (date_dif == min(date_dif));  % pick the windspeed from the neareast field trip date when we have the windspeed data
end

%% data quality control
% the wind speed on Jan 20th, 2019 (12.94m/s) looks suspicious, showing
% large difference from its neighboring dates, therefore I replace that
% value with the mean value of the wind speed on the neighboring dates

Kusu.wind_speed (29) = mean([9.0564 9.1974 9.45 9.0173]);

%%
% F = k660 K0 (CO2_water - CO2_air)
% pCO2_air from NOAA ESRL Carbon Cycle Coop- erative Global Air Sampling Network (Dlugokencky et al., 2018)
% k660 from Wanninkhof et al. 2009
% K0 from Weiss 1974, lnK0 = A1+A2*(100/T) + A3 * ln(T/100) + S* (B1 + B2/(T/100) + B3*(T/100)^2)
% T is absolute temperature

% function format: K0 = sol_calc (T, Sal, mode)         k660 (cm/hr) = k_calc (T,
% wind speed (unit:m/s));

T = fulldata_co2.T;  %temperature, unit:Celcius
salinity = fulldata_co2.sal;
u = fulldata_co2.wind_speed;  % wind speed, unit: m/s, or we can use the value from each day
pCO2_w = fulldata_co2.co2 .* 10^-6;
pCO2_a = 400*(10^-6);   %unit:atm, pCO2 in the air, Wit et al. 2018

K0 = sol_calc (T, salinity, 2);  %calculate the solubility constant, input temperature, salinity and unit (1 means mol/kg.atm; 2 means mol/L.atm)
k660 = k_calc (T, u);    %calculate the gas transfer velocity using schmidt number, input temperature and wind speed

% calculate the CO2 flux for each day whenever we have pCO2 data
fulldata_co2.flux = k660 .* K0 .* (pCO2_w - pCO2_a) .* 10 .* 24;    %unit: mol/(day*m2)

%% plot
figure('color','w');
plot(fulldata_co2.date, fulldata_co2.flux, 'ko','markerfacecolor','k')
ylabel('CO2 flux (mol/day/m2)')
xlabel('Time')

figure('color','w');
plot(fulldata_co2.date, fulldata_co2.wind_speed, 'ko','markerfacecolor','k')
ylabel('wind speed(m/s)')
xlabel('Time')


figure('color','w');
plot(fulldata_co2.date, fulldata_co2.co2, 'ko','markerfacecolor','k')
ylabel('pCO2(uatm)')
xlabel('Time')

%% calculate the annual CO2 flux
startdate = datetime('1-jan-2019');
enddate = datetime('9-jan-2020');

subjoey = fulldata_co2 (fulldata_co2.date > startdate & fulldata_co2.date < enddate & ~isnan(fulldata_co2.flux),:);
subjoey.datenum = round(datenum( subjoey.date - startdate));
figure;plot(subjoey.date, subjoey.flux, 'k-')
co2_flux = trapz(subjoey.datenum, subjoey.flux);  % mol/m2/yr
co2_flux_areatotal = co2_flux * 12 * 127674 * (1000)^2 / (10)^12; % Tg/yr


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



%%
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\pCO2\fulldata_co2.mat');
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat','Kusu');

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

function [k660] = k_calc(T, u)   %T is in degree Celcius   (Wanninkhof, 1992; Jahne et al. 1987b)
    A = 2073.1;
    B = 125.62;
    C = 3.6276;
    D = 0.043219;
    Sc = A - B*T + C*T.^2 - D*T.^3;  %Schmidt number dependent on temperature
    k660 = 0.39*u.^2 .* (Sc/660).^(-0.5);
end