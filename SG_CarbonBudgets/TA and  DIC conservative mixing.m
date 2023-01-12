% set up the endmembers
clear
clear global model
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\endmember.mat')
global model

% load the Kusu Hantu data
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat','Kusu','Hantu');
load ('F:\NTU\Research\Kusu Hantu Biogeochem\Seabird Sal\processed data for pH correction\daily mean salinity\dailymean_sal.mat','dm_sal');




%%
figure('color','w')
plot(model.sal, DIC_cons_calc(model.sal), 'k-');hold on   % plot the conservative mixing line
plot(Kusu.S_Valeport, Kusu.DIC,'ko','markerfacecolor','k');hold on
plot(Hantu.S_Valeport, Hantu.DIC, 'ko','markerfacecolor','k');hold on
plot(Kusu.S_Valeport(strcmp(Kusu.season,'SW monsoon')==1), Kusu.DIC(strcmp(Kusu.season,'SW monsoon')==1), 'go','markerfacecolor','g');hold on
plot(Hantu.S_Valeport(strcmp(Hantu.season,'SW monsoon')==1), Hantu.DIC(strcmp(Hantu.season,'SW monsoon')==1), 'go','markerfacecolor','g');hold on
plot(Kusu.S_Valeport(strcmp(Kusu.season,'NE monsoon')==1), Kusu.DIC(strcmp(Kusu.season,'NE monsoon')==1), 'bo','markerfacecolor','b');hold on
plot(Hantu.S_Valeport(strcmp(Hantu.season,'NE monsoon')==1), Hantu.DIC(strcmp(Hantu.season,'NE monsoon')==1), 'bo','markerfacecolor','b');hold on
 
 
 
xlim([29 34])
ylabel('DIC (umol/kg)')
xlabel('salinity')
 
 
%% TA
% conservative mixing model
figure('color','w')
plot(model.sal, TA_cons_calc(model.sal), 'k-');hold on   % plot the conservative mixing line
plot(Kusu.S_Valeport, Kusu.TA,'ko','markerfacecolor','k');hold on
plot(Hantu.S_Valeport, Hantu.TA, 'ko','markerfacecolor','k');hold on
plot(Kusu.S_Valeport(strcmp(Kusu.season,'SW monsoon')==1), Kusu.TA(strcmp(Kusu.season,'SW monsoon')==1), 'go','markerfacecolor','g');hold on
plot(Hantu.S_Valeport(strcmp(Hantu.season,'SW monsoon')==1), Hantu.TA(strcmp(Hantu.season,'SW monsoon')==1), 'go','markerfacecolor','g');hold on
plot(Kusu.S_Valeport(strcmp(Kusu.season,'NE monsoon')==1), Kusu.TA(strcmp(Kusu.season,'NE monsoon')==1), 'bo','markerfacecolor','b');hold on
plot(Hantu.S_Valeport(strcmp(Hantu.season,'NE monsoon')==1), Hantu.TA(strcmp(Hantu.season,'NE monsoon')==1), 'bo','markerfacecolor','b');hold on
 
xlim([29 34])
ylabel('TA(umol/kg)')
xlabel('salinity')

%% calculate the conservative mixing DIC and TA
Kusu.DIC_cons = DIC_cons_calc(Kusu.salinity);
Hantu.DIC_cons = DIC_cons_calc(Hantu.salinity);
Kusu.TA_cons = TA_cons_calc(Kusu.salinity);
Hantu.TA_cons = TA_cons_calc(Hantu.salinity);





%% time series of deviation from conservative mixing
% %DIC 
Kusu.DIC_dev = Kusu.DIC - Kusu.DIC_cons;
Hantu.DIC_cons = DIC_cons_calc(Hantu.salinity);
Hantu.DIC_dev = Hantu.DIC - Hantu.DIC_cons;
%
figure('color','w')
plot([datetime('01-may-2017') datetime('31-Dec-2020')], [0 0],'k-','linewidth',2);hold on
plot(Kusu.Date,Kusu.DIC_dev,'bo','markerfacecolor','b','markersize',8);hold on
plot(Hantu.Date,Hantu.DIC_dev,'bo','markerfacecolor','b','markersize',8);

xlim([datetime('01-Oct-2017') datetime('31-Dec-2020')])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
YTICK = [-40 -20 0 20 40]
set(gca,'xtick',datetime(XTICK), 'xticklabel',datestr(XTICK,'mmm-yy'),'YTICK',YTICK)
ylabel('Deviation of DIC from conservative mixing (umol kg^-^1)');
xlabel('Time')
set(gca,'fontsize',14)


%% time series of deviation from conservative mixing
% TA
Kusu.TA_dev = Kusu.TA - Kusu.TA_cons;
Hantu.TA_dev = Hantu.TA -  Hantu.TA_cons;
 
figure('color','w')
plot([datetime('01-may-2017') datetime('31-Dec-2020')], [0 0],'k-','linewidth',2);hold on
plot(Kusu.Date,Kusu.TA_dev,'bo','markerfacecolor','b');hold on
plot(Hantu.Date,Hantu.TA_dev,'bo','markerfacecolor','b');

xlim([datetime('01-Oct-2017') datetime('31-dec-2020')])
ylim([-44 40])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
YTICK = [-40 -20 0 20 40];
set(gca,'xtick',datetime(XTICK), 'xticklabel',datestr(XTICK,'mmm-yy'),'ytick',YTICK);
ylabel('Deviation of TA from conservative mixing (umol kg^-^1)');
xlabel('Time')
set(gca,'fontsize',14)
 
%% scatter plot of deviation of DIC vs TA 

figure('color','w')

h1=plot(Kusu.DIC_dev, Kusu.TA_dev , 'o', 'markeredgecolor',[102 179 255]/255 ,'markerfacecolor',[102 179 255]/255);hold on
%h2=plot(Kusu.DIC_dev (strcmp(Kusu.season,'NE monsoon')==1), Kusu.TA_dev (strcmp(Kusu.season,'NE monsoon')==1), 'o','markerfacecolor',[92 167 186]/255 , 'markeredgecolor');hold on
h3=plot(Kusu.DIC_dev (strcmp(Kusu.season,'SW monsoon')==1), Kusu.TA_dev (strcmp(Kusu.season,'SW monsoon')==1), 'o','markerfacecolor',[107 194 53]/255, 'markeredgecolor', [107 194 53]/255);hold on


plot(Hantu.DIC_dev, Hantu.TA_dev, 'o','markeredgecolor',[102 179 255]/255 ,'markerfacecolor',[102 179 255]/255);hold on
%plot(Hantu.DIC_dev (strcmp(Hantu.season,'NE monsoon')==1), Hantu.TA_dev (strcmp(Hantu.season,'NE monsoon')==1), 'bo','markerfacecolor','b');hold on
plot(Hantu.DIC_dev (strcmp(Hantu.season,'SW monsoon')==1), Hantu.TA_dev (strcmp(Hantu.season,'SW monsoon')==1), 'o','markerfacecolor',[107 194 53]/255, 'markeredgecolor', [107 194 53]/255);hold on

%plot the reference lines for biogeochem processes
plot([0 30], [0 -1/40*30], 'k-');hold on   %remineralization, use the C/N ratio of 40 from Baum 2008
plot([0 -30],[0 16/106*30], 'k-');hold on   %photosynthesis, use redfield ratio
plot([0 15],[0 30],'k-');hold on
plot([0 -15],[0 -30],'k-');

YTICK = [-40 -20 0 20 40];
XTICK = [-40 -20 0 20 40];
ylim([-50 50])
xlim([-50 50])
xlabel ('DIC deviation')
ylabel  ('TA deviation')
set(gca,'fontsize',14,'ytick',YTICK, 'XTICK', XTICK)
legend([h1  h3],'other seasons','SW monsoon season')


%% time series of DIC

%% conservative mixing line for DIC
% 
% import seabird CTD data
load('F:\NTU\Research\Kusu Hantu Biogeochem\Seabird Sal\processed data for pH correction\daily mean salinity\dailymean_sal.mat')
dm_sal.salinity(807:837) = NaN;  % these data are believed to be wrong due to large discrepancy with the Valeport sensor (04Nov2018 to 04Dec2018)

% calculate conservative mixing DIC and TA using seabird data
f_riv = 1 - Kusu.salinity ./ model.sal_mar;
f_mar = Kusu.salinity ./ model.sal_mar;
DIC_cons = DIC_cons_calc(dm_sal.salinity);
TA_cons = TA_cons_calc(dm_sal.salinity);


%%
figure('color','w');
col_cons = [97 227 85]/255;
col_meas = [102 179 255]/255;
plot(dm_sal.DT, DIC_cons, 'k-');hold on
plot(Kusu.Date, Kusu.DIC, 'o', 'markeredgecolor',col_meas,'markerfacecolor', col_meas,'markersize',8);hold on
plot(Hantu.Date, Hantu.DIC, 'o', 'markeredgecolor',col_meas,'markerfacecolor', col_meas,'markersize',8)
xlim([datetime('01-Oct-2017') datetime('31-Dec-2020')])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
set(gca,'xtick',datetime(XTICK), 'xticklabel',datestr(XTICK,'mmm-yy'))
%ylim([-2.2 0.8])
%ylabel('d13C-DIC (‰)')
xlabel('Time')
set(gca,'fontsize',14)



%% save data
save ('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\DOC mineralizatoin\DIC_TA_deviations.mat');
save('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat','Kusu','Hantu')


%% nested functions
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



function d13c_cons = d13c_cons_calc (salinity)
    global model
    f_river = 1 - salinity ./ model.sal_mar;
    f_marine = salinity ./ model.sal_mar;
    A = 1./ (1 ./ model.R_riv + 1);
    B = 1./ (1 ./ model.R_mar + 1);
    C = 1./ (model.R_riv + 1);
    D = 1./ (model.R_mar + 1);
    c13ratio_cons = (A .* f_river .* model.dic_riv + B .* f_marine .* model.dic_mar) ./ (C .* f_river .* model.dic_riv + D .* f_marine .* model.dic_mar) ;
    d13c_cons = (c13ratio_cons - model.R_vpdb) ./ model.R_vpdb .* 1000;
    
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


