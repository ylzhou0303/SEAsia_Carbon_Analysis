%% Use a new approach that compiles the equations of DIC_dev, TA_dev and the d13C-DIC as an equation group
% to calculate the molar contribution of remineralization, dissolution and
% CO2 outgassing

% import data
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat','Kusu','Hantu')
clear global model
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\endmember.mat')
global model

%% calculate the amount of tDOC
model.R_peat = model.R_vpdb * (model.d13c_peat ./ 1000 + 1);
model.R_DOC_mar = model.d13cDOC_mar ./1000 .* model.R_vpdb + model.R_vpdb;   %marine endmember of d13c-DOC

Kusu.R_DOC_meas = Kusu.d13cDOC ./1000 .* model.R_vpdb + model.R_vpdb;  
Kusu.tDOC = (Kusu.R_DOC_meas - model.R_DOC_mar) .* Kusu.DOC ./ (model.R_peat - model.R_DOC_mar);
Kusu.tDOC_perc = Kusu.tDOC ./ Kusu.DOC *100;
Kusu.tDOC(Kusu.tDOC<0) = NaN; %negative result means no tDOC existing, but it's not good to say it's zero

Hantu.R_DOC_meas = Hantu.d13cDOC ./1000 .* model.R_vpdb + model.R_vpdb;  
Hantu.tDOC = (Hantu.R_DOC_meas - model.R_DOC_mar) .* Hantu.DOC ./ (model.R_peat - model.R_DOC_mar);
Hantu.tDOC_perc = Hantu.tDOC ./ Hantu.DOC *100;
Hantu.tDOC(Hantu.tDOC<0) = NaN;

%convert the tDOC and mDOC to umol/kg because I need to add tDOC with tDIC,
%and the latter is umol/kg
Kusu.mDOC = Kusu.DOC - Kusu.tDOC;
Hantu.mDOC = Hantu.DOC - Hantu.tDOC;
Kusu.tDOC = unitqlo(Kusu.tDOC, Kusu.salinity, Kusu.temperature, 1);
Hantu.tDOC = unitqlo(Hantu.tDOC, Hantu.salinity, Hantu.temperature, 1);
Kusu.mDOC = unitqlo(Kusu.mDOC, Kusu.salinity, Kusu.temperature, 1);
Hantu.mDOC = unitqlo(Hantu.mDOC, Hantu.salinity, Hantu.temperature, 1);


%% make a plot for the %tDOC
figure('color','w')
plot(Kusu.Date, Kusu.tDOC ./ Kusu.DOC *100, 'bo', 'markerfacecolor', 'b');hold on
plot(Hantu.Date, Hantu.tDOC ./ Hantu.DOC *100,  'bo', 'markerfacecolor', 'b');

ylabel('%tDOC')
xlim([datetime('01-Oct-2017') datetime('31-Dec-2020')])
YTICK = [0 10 20 30 40 50]
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
set(gca,'xtick',XTICK,'xticklabel',datestr(XTICK,'mmm-yy'),'ytick',YTICK,'fontsize',14)


%% isotopic mass balance calculation
e = 23.644 - 9701.5 ./ (Kusu.temperature + 273.15); %Rau et al. 1996
a = e ./ 1000 + 1; % calculate the fractionation factor a for outgassing


Kusu.REM = zeros(size(Kusu,1),1);
Kusu.REM (Kusu.REM == 0) = NaN;
A = 1 ./ Kusu.DIC_cons .* (model.d13c_tDIC - Kusu.d13cDIC_cons); 
B = 1 ./ Kusu.DIC_cons .* 1000 .* (a - 1);
C = 1 ./ Kusu.DIC_cons .* (model.d13c_carb - Kusu.d13cDIC_cons);
        
Kusu.REM = ((Kusu.d13cDIC - Kusu.d13cDIC_cons) - B.*Kusu.DIC_dev + 0.5 .* Kusu.TA_dev .* (C - B)) ./ (A - 1.0125 .* B + 0.0125 .* C);
Kusu.DIS = 0.5 .* Kusu.TA_dev + 0.0125*Kusu.REM;
Kusu.CO2 = Kusu.REM + Kusu.DIS - Kusu.DIC_dev;

%% for Hantu
e = 23.644 - 9701.5 ./ (Hantu.temperature + 273.15); %Rau et al. 1996
a = e ./ 1000 + 1;
 

Hantu.REM = zeros(size(Hantu,1),1);
Hantu.REM (Hantu.REM == 0) = NaN;
A = 1 ./ Hantu.DIC_cons .* (model.d13c_tDIC - Hantu.d13cDIC_cons); 
B = 1 ./ Hantu.DIC_cons .* 1000 .* (a - 1);
C = 1 ./ Hantu.DIC_cons .* (model.d13c_carb - Hantu.d13cDIC_cons);
        
Hantu.REM = ((Hantu.d13cDIC - Hantu.d13cDIC_cons) - B.*Hantu.DIC_dev + 0.5 .* Hantu.TA_dev .* (C - B)) ./ (A - 1.0125 .* B + 0.0125 .* C);
Hantu.DIS = 0.5 .* Hantu.TA_dev + 0.0125*Hantu.REM;
Hantu.CO2 = Hantu.REM + Hantu.DIS - Hantu.DIC_dev;

%% monthly average
% compile data from Kusu and Hantu into one single dataset
joey = table();
joey.date = vertcat(Kusu.Date, Hantu.Date);
joey.temperature = vertcat(Kusu.temperature, Hantu.temperature);
joey.location = vertcat(Kusu.location, Hantu.location);
joey.season = vertcat(Kusu.season, Hantu.season);
joey.salinity = vertcat(Kusu.salinity, Hantu.salinity);
joey.d13cDOC = vertcat(Kusu.d13cDOC, Hantu.d13cDOC);
joey.DOC = vertcat(Kusu.DOC, Hantu.DOC); % unit:umol/L
joey.mDOC = vertcat(Kusu.mDOC, Hantu.mDOC); % unit: umol/kg
joey.tDOC = vertcat(Kusu.tDOC, Hantu.tDOC); % unit: umol/kg
joey.REM = vertcat(Kusu.REM, Hantu.REM);% unit:umol/kg
joey.CO2 = vertcat(Kusu.CO2, Hantu.CO2);
joey.DIS = vertcat(Kusu.DIS, Hantu.DIS);
joey.tDOC_perc = vertcat(Kusu.tDOC_perc, Hantu.tDOC_perc);

joey.REM (joey.REM < 0) = 0;
joey.CO2 (joey.CO2 < 0) = 0;
joey.CO2 (joey.REM < 0) = 0;
joey.CO2 (joey.CO2 > joey.REM) = joey.REM(joey.CO2 > joey.REM);
joey.REM_remain = joey.REM - joey.CO2;
% negative number means no tDOC remineralization happening, primary
% production happens instead. here we only discuss remieralization, so set
% the negative numbers as 0
% similarly, negative CO2 means no CO2 outgassing happening. Instead, CO2
% uptake happends. Here we only show the intensity of CO2 outgassing, so
% set the negative numbers as 0.

%% calculate monthly mean
mthave = table();
mth = datetime('15-jan-2019');
i=1;
while datenum(mth - joey.date(end))<30
    mthave.date(i) = mth;
    p = find (strcmp(cellstr(datestr(joey.date,'mmm-yyyy')),datestr(mth,'mmm-yyyy')) == 1);
    %look for the data for that month
    
    mthave.tDOC(i) = mean(joey.tDOC(p),'omitnan');
    mthave.REM(i) = mean(joey.REM(p),'omitnan');
    mthave.CO2(i) = mean(joey.CO2(p),'omitnan');
    mthave.REM_remain(i) = mean(joey.REM_remain(p),'omitnan');
    
    mth = mth + datenum(30);
    i = i+1;
end

mthave.REM_perc = mthave.REM ./ (mthave.REM + mthave.tDOC);

for i = 1:size(mthave,1)
    tempchar = datestr(mthave.date(i));
    tempchar = ['15',tempchar(3:end)];
    mthave.date(i) = cellstr(tempchar);
end


%% GREAT!Then make the barplots

figure('color','w');
y_bar = horzcat(mthave.tDOC, mthave.REM_remain, mthave.CO2);
p = bar(mthave.date, y_bar, 0.65, 'stacked', 'facecolor','flat');hold on
p(1).FaceColor = [69 137 148]/255;
p(2).FaceColor = [178 200 187]/255;
p(3).FaceColor = [102 179 255]/255;

ylabel('terrigenous carbon budget (umol/L)')

XTICKS = datetime({'15-Jan-2019','15-Apr-2019','15-Jul-2019','15-Oct-2019','15-Jan-2020','15-Apr-2020','15-Jul-2020','15-Oct-2020'});
set(gca,'xtick',XTICKS,'xticklabel',datestr(XTICKS,'mmm-yyyy'))
legend('tDOC remaining','tDOC remineralized as DIC','tDOC remineralized and degassed');
% ylabel('percentage of tDOC (%)')
ylabel('carbon budget (umol/L)')
xlabel('Time')
xlim([datetime('15-Dec-2018') datetime('31-dec-2020')])
set(gca,'fontsize',14);


%% extrapolation for the riverine endmember
joey.total_tDOC  = joey.tDOC + joey.REM;

X = joey.salinity;
X = horzcat(ones(size(X,1),1), X);
Y = joey.total_tDOC;
[B,BINT,R,RINT,STATS] = regress(Y,X);

YTICK = [0 150 300 450 600 750 900];
figure('color','w');
plot([0:35],B(1)+B(2).*[0:35],'k-');hold on
plot(X(:,2),Y,'bo','markerfacecolor','b')
set(gca,'ytick',YTICK);

%% extrapolation for the riverine endmember (with only SW Monsoon and inter-monsoon)
joey2 = joey(strcmp(joey.season,'NE monsoon') == 0, :);
joey2.total_tDOC  = joey2.tDOC + joey2.REM;

X = joey2.salinity;
X = horzcat(ones(size(X,1),1), X);
Y = vertcat(joey2.total_tDOC);
Y = unitqlo(Y, joey2.salinity, joey2.temperature, 2);  %convert it to umol/L
[B,BINT,R,RINT,STATS] = regress(Y,X);


%% close-up
figure('color','w');
plot([0:35],B(1)+B(2).*[0:35],'k-');hold on
plot(X(:,2),Y,'bo','markerfacecolor','b');
errorbar(0, B(1), 104, 'bs','markerfacecolor', 'b');
xlim([28 35])
ylim([-20 160])
%YTICK = [0 150 300 450 600 750 900 1050];
YTICK = [0 30 60 90 120 150];
set(gca,'ytick',YTICK, 'fontsize', 14);
xL=xlim;yL=ylim;
plot(xL,[yL(2),yL(2)],'k',[xL(2),xL(2)],[yL(1),yL(2)],'k')
box off
axis([xL yL])


%% full salinity
figure('color','w');
plot([0:35],891.6+(-26.5).*[0:35],'k-');hold on
plot(X(:,2),Y,'b.','markerfacecolor','b','markersize',10);
errorbar(0, 891.63, 104, 'bs','markerfacecolor', 'b');
errorbar(-2, 892, 159, 'rs', 'markerfacecolor', 'r');
xlim([-4 35])
ylim([-20 1100])
YTICK = [0 150 300 450 600 750 900 1050];
set(gca,'ytick',YTICK, 'fontsize', 14);
xL=xlim;yL=ylim;
plot(xL,[yL(2),yL(2)],'k',[xL(2),xL(2)],[yL(1),yL(2)],'k')
box off
axis([xL yL])


%% use the results to reconstruct the DIC, TA and d13C-DIC
Kusu.DIC_mod = Kusu.DIC_cons + Kusu.REM + Kusu.DIS + (-1) * Kusu.CO2;
Kusu.TA_mod = Kusu.TA_cons + (-0.025) * Kusu.REM + 2*Kusu.DIS;
Kusu.d13cDIC_mod = Kusu.d13cDIC_cons + Kusu.REM ./ Kusu.DIC_cons .* (model.d13c_tDIC - Kusu.d13cDIC_cons) + Kusu.DIS ./ Kusu.DIC_cons .* (model.d13c_carb - Kusu.d13cDIC_cons) + (-Kusu.CO2 ./ Kusu.DIC_cons) .* 1000 .* (a - 1);

figure;
plot(Kusu.Date, Kusu.d13cDIC, 'ko');hold on
plot(Kusu.Date, Kusu.d13cDIC_mod, 'ro')
figure;
plot(Kusu.d13cDIC, Kusu.d13cDIC_mod, 'ko'); hold on
plot([-2 2],[-2 2],'k-')

%%
figure;
plot(Kusu.Date,Kusu.d13cDIC_cons,'go');hold on
plot(Kusu.Date,Kusu.d13cDIC_cons + Kusu.REM ./ Kusu.DIC_cons .* (model.d13c_tDIC - Kusu.d13cDIC_cons), 'ro');hold on
%plot(Kusu.Date,Kusu.d13cDIC_cons + Kusu.DIS ./ Kusu.DIC_cons .* (model.d13c_carb - Kusu.d13cDIC_cons), 'wo');hold on
plot(Kusu.Date,Kusu.d13cDIC_cons + (-Kusu.CO2 ./ Kusu.DIC_cons) .* 1000 .* (a - 1),'bo');hold on
plot(Kusu.Date,Kusu.d13cDIC, 'ko')

legend('cons','REM','CO2','measured')

%%
save('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\DOC mineralizatoin\tDOC_budget_results.mat')
save('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat','Kusu','Hantu');
%% nested functions
function [a_DIC_g] = a_DIC_g_calc(T,pH)
    T (T>25) = 25;  %because this equation is only applicable up to 25C
    e_DIC_g = 0.0144 * T .* f_CO3_calc(T,pH) - 0.107 * T + 10.53;
    a_DIC_g = e_DIC_g / 1000 + 1;
end

function [f_CO3] = f_CO3_calc(T,pH)
    T (T>25) = 25;  %because the equation above is only applicable up to 25C
    T = T + 273.15;
    Ka1 = 10.^(-(3404.71 ./ T + 0.032786 .* T - 14.8435));
    Ka2 = 10.^(-(2902.39 ./ T + 0.02379 .* T - 6.4980));
    f_CO3 = Ka1 .* Ka2 ./ ((10.^ (-pH)).^2 + 10.^ (-pH) .* Ka1 + Ka1 .* Ka2);
end


function R = R_calc (d13c)
    global model
    R = (d13c ./ 1000 + 1) .* model.R_vpdb;
end

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