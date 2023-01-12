%% calculate the respective contribution to the pH (pH budget)
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat','Kusu','Hantu');
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\DOC mineralizatoin\tDOC_budget_results.mat','joey')

%% calculate the pH due to conservative mixing
% calculate the DIC and TA due to calcium carbonate dissolution/production
Kusu.silicate (isnan(Kusu.silicate)) = 5;
Kusu.phosphate (isnan(Kusu.phosphate)) = 0.13;


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
A=CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);
Kusu.pH_cons = A(:,3);
Kusu.Hfree_cons = A(:,13);





%% calculate the pH due to calcium carbonate dissolution/production
% calculate the DIC and TA due to calcium carbonate dissolution/production
Kusu.silicate (isnan(Kusu.silicate)) = 5;
Kusu.phosphate (isnan(Kusu.phosphate)) = 0.13;


par1type =    1; % The first parameter supplied is of type "1", which is "alkalinity"
par1     =   Kusu.TA_cons + 2*Kusu.DIS; % value of the first parameter
par2type =    2; % The 2nd parameter supplied is of type "2", which is "DIC"
par2     =   Kusu.DIC_cons + Kusu.DIS; % value of the second parameter, which is a long vector of different DIC's!
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
Kusu.pH_DIS = A(:,3);
Kusu.Hfree_DIS = A(:,13);

%% calculate the pH due to remineralization of natural tDOC flux
Kusu.silicate (isnan(Kusu.silicate)) = 5;
Kusu.phosphate (isnan(Kusu.phosphate)) = 0.13;

par1type =   1; % The first parameter supplied is of type "1", which is "alkalinity"
par1     =   Kusu.TA_cons + 2 * Kusu.DIS + (-0.025).* Kusu.REM .* 0.65 .* (Kusu.REM>0) + (-0.0625) .* Kusu.REM .* (Kusu.REM<0); % value of the first parameter
% only apply the 0.65 (natural fraction) to the dates when we have net
% remineralization, we consider the net production of OM is all natural
% process
par2type =    2; % The 2nd parameter supplied is of type "2", which is "DIC"
par2     =   Kusu.DIC_cons + Kusu.DIS + (Kusu.REM + (-1)*Kusu.CO2) .* 0.65 .* (Kusu.REM>0) + (Kusu.REM + (-1)*Kusu.CO2) .* (Kusu.REM<0); % value of the second parameter, which is a long vector of different DIC's!
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
Kusu.pH_DIS_NATREM = A(:,3);
Kusu.Hfree_DIS_NATREM = A(:,13);




%% calculate the pH due to calcium carbonate dissolution/production + natural + anthropogenic remineralization
% calculate the DIC and TA due to calcium carbonate dissolution/production
Kusu.silicate (isnan(Kusu.silicate)) = 5;
Kusu.phosphate (isnan(Kusu.phosphate)) = 0.13;


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
Kusu.pH_ANT = A(:,3);
Kusu.Hfree_ANT = A(:,13);



%%
figure;
plot(Kusu.Date, Kusu.pH_cons,'go');hold on
plot(Kusu.Date, Kusu.pH_DIS, 'bo');hold on
plot(Kusu.Date, Kusu.pH_DIS_NATREM,'ro');hold on
plot(Kusu.Date, Kusu.pH_calc,'ko');

legend('cons','calcium carbonate','natural tDOC','natural + human')


%%
figure;
plot(Kusu.Date, Kusu.Hfree_cons,'go');hold on
plot(Kusu.Date, Kusu.Hfree_DIS, 'bo');hold on
plot(Kusu.Date, Kusu.Hfree_DIS_NATREM,'ro');hold on
plot(Kusu.Date, Kusu.Hfree_ANT,'ko');

legend('cons','calcium carbonate','natural tDOC','natural + human')


%% calculate the pH budget and the monthly average
Kusu.H1 =  Kusu.Hfree_DIS - Kusu.Hfree_cons;
Kusu.H2 = Kusu.Hfree_DIS_NATREM - Kusu.Hfree_DIS; % the contribution of natural tDOC flux to the pH reduction
Kusu.H3 = Kusu.Hfree_ANT - Kusu.Hfree_DIS_NATREM;

H_mthave = table();
mth = datetime('20-jan-2018');
i=1;
while datenum(mth - Kusu.Date(end))<30
    H_mthave.date(i) = mth;
    p = find (strcmp(cellstr(datestr(Kusu.Date,'mmm-yyyy')),datestr(mth,'mmm-yyyy')) == 1);
    %look for the data for that month
    
    H_mthave.H1(i) = mean(Kusu.H1(p),'omitnan');
    H_mthave.H2(i) = mean(Kusu.H2(p),'omitnan');
    H_mthave.H3(i) = mean(Kusu.H3(p),'omitnan');
    
    mth = mth + datenum(30);
    i = i+1;
end

%%
H_mthave.ANT_perc = H_mthave.H3./(H_mthave.H1 + H_mthave.H2 + H_mthave.H3);

figure;
plot(H_mthave.date, H_mthave.H3,'ko')
%%
for i = 1:size(ph_mthave,1)
    tempchar = datestr(ph_mthave.date(i));
    tempchar = ['15',tempchar(3:end)];
    ph_mthave.date(i) = cellstr(tempchar);
end


%% make bar plot
H_mthave_sw = H_mthave([5:9 17:21 29:33],:);  %only show the SW Monsoon season
y_bar = horzcat(H_mthave_sw.H1, H_mthave_sw.H2, H_mthave_sw.H3);
y_bar_pos = y_bar;
y_bar_pos (y_bar_pos < 0) = 0;
y_bar_neg = y_bar;
y_bar_neg (y_bar_neg > 0) = 0;

figure('color','w');
% we just only show the increase in H concentration
p2 = bar(H_mthave_sw.date, y_bar_pos .* 1000, 0.8, 'stacked', 'facecolor','flat');  %unit: nmol/kg
p2(1).FaceColor = [250	250	250]/255;
p2(2).FaceColor = [198	229	243]/255;
p2(3).FaceColor = [83	157	219]/255;


xlim([datetime('01-Oct-2017') datetime('31-Dec-2020')])
%ylim([-0.13 0])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
set(gca,'xtick',datetime(XTICK), 'xticklabel',datestr(XTICK,'mmm-yy'))
set(gca,'fontsize',14)

%%
save('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\DOC mineralizatoin\phbudget.mat');
save('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat','Kusu','Hantu');


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