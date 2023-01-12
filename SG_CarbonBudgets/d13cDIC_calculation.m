%% import data
clear global model
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\endmember.mat')
global model

load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat','Kusu','Hantu')

% import seabird CTD data
load('F:\NTU\Research\Kusu Hantu Biogeochem\Seabird Sal\processed data for pH correction\daily mean salinity\dailymean_sal.mat')
dm_sal.salinity(807:837) = NaN;  % these data are believed to be wrong due to large discrepancy with the Valeport sensor (04Nov2018 to 04Dec2018)

%% calculate the conservative mixing d13C-DIC, using the Kusu Hantu data
f_riv = 1 - Kusu.salinity ./ model.sal_mar;
f_mar = Kusu.salinity ./ model.sal_mar;
Kusu.d13cDIC_cons = d13c_mix (model.dic_riv, model.d13c_riv, f_riv, model.dic_mar, model.d13c_mar, f_mar);

f_riv = 1 - Hantu.salinity ./ model.sal_mar;
f_mar = Hantu.salinity ./ model.sal_mar;
Hantu.d13cDIC_cons = d13c_mix (model.dic_riv, model.d13c_riv, f_riv, model.dic_mar, model.d13c_mar, f_mar);


%%  calculate the conservative mixing d13C-DIC, using the seabird CTD data
f_riv = 1 - dm_sal.salinity ./ model.sal_mar;
f_mar = dm_sal.salinity ./ model.sal_mar;
d13cDIC_cons_CTD = d13c_mix (model.dic_riv, model.d13c_riv, f_riv, model.dic_mar, model.d13c_mar, f_mar);

%% plot
Hantu.d13cDIC([49:52]) = NaN;

figure('color','w');
col_cons = [97 227 85]/255;
col_meas = [102 179 255]/255;
plot(d13cDIC.DT, d13cDIC.cons, 'o-','markeredgecolor',col_cons,'markerfacecolor',col_cons,'markersize',5);hold on
plot(Kusu.Date, Kusu.d13cDIC, 'o', 'markeredgecolor',col_meas,'markerfacecolor', col_meas,'markersize',8);hold on
plot(Hantu.Date, Hantu.d13cDIC, 'o', 'markeredgecolor',col_meas,'markerfacecolor', col_meas,'markersize',8)
xlim([datetime('01-Oct-2017') datetime('31-Dec-2020')])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
set(gca,'xtick',datetime(XTICK), 'xticklabel',datestr(XTICK,'mmm-yy'))
ylim([-2.2 0.8])
ylabel('d13C-DIC (‰)')
xlabel('Time')
set(gca,'fontsize',14)

%% save data
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\DOC mineralizatoin\d13cDIC.mat');
save('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat','Kusu','Hantu');
%% nested functions
function [d13c_mixed] = d13c_mix (pool_A, d13c_A, f_A, pool_B, d13c_B, f_B)   
%calculate the resulting d13c if two pools of DIC is mixed together
   
    global model
    R_A = (d13c_A ./ 1000 + 1) .* model.R_vpdb;
    R_B = (d13c_B ./ 1000 + 1) .* model.R_vpdb;
    
    C13 = pool_A .* R_A ./(1+R_A) .* f_A + pool_B .* R_B ./(1+R_B) .* f_B;
    C12 = pool_A .* R_A ./(1+R_A) ./ R_A .* f_A + pool_B .* R_B ./(1+R_B) ./ R_B .* f_B;
    
    R_mix = C13 ./ C12;
    d13c_mixed = (R_mix ./ model.R_vpdb - 1) .* 1000;
end