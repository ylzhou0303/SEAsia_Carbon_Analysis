%% uncertainty analysis based on monte carlo simulation
% IN THIS SCRIPT, I CALCULATE THE MONTHLY AVERAGE FOR ALL PARAMETERS
% INVOLVED FIRST FOR MONTE CARLO SIMULATION

save('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat','Kusu','Hantu');

clear global model
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\endmember.mat')
global model

MC = struct();
%%
rep = 10000;  % iteration: 1000 times for Monte Carlo Simulation

m = size(Kusu,1) + size(Hantu, 1);
%generate a matrix of uncertainty for the input variables
% (1) for marine and riverine endmebers
MC.dic_mar = randn(1,rep) * model_std.dic_mar_std  + model.dic_mar; MC.dic_mar = repmat (MC.dic_mar, m, 1);
MC.ta_mar = randn(1,rep) * model_std.ta_mar_std  + model.ta_mar; MC.ta_mar = repmat (MC.ta_mar, m, 1);
MC.sal_mar = randn(1,rep) * model_std.sal_mar_std  + model.sal_mar; MC.sal_mar = repmat (MC.sal_mar, m, 1);
MC.d13c_mar = randn(1,rep) * model_std.d13cdic_mar_std  + model.d13c_mar; MC.d13c_mar = repmat (MC.d13c_mar, m, 1);
MC.d13cdoc_mar = randn(1,rep) * model_std.d13cDOC_mar_std  + model.d13cDOC_mar; MC.d13cdoc_mar = repmat (MC.d13cdoc_mar, m, 1);

MC.dic_riv = randn(1,rep) * 34 + model.dic_riv; MC.dic_riv = repmat (MC.dic_riv, m, 1);
MC.ta_riv = randn(1,rep) * 34 + model.ta_riv; MC.ta_riv = repmat (MC.ta_riv, m, 1);
MC.d13c_peat = randn(1,rep) * 1  + model.d13c_peat;  %1 permil is 1SD 
MC.d13c_peat = repmat (MC.d13c_peat, m, 1);
MC.d13c_tDIC = randn(1,rep) * 1 + model.d13c_tDIC; MC.d13c_tDIC = repmat (MC.d13c_tDIC, m, 1);
MC.d13c_riv = randn(1,rep) * 1 + model.d13c_riv; %assume a SD of 1 permil for the riverine endmember
MC.d13c_riv = repmat (MC.d13c_riv, m, 1);
 
% %% (2) for the actual measurements (use the standard deviation of the triplicate measurements)
% for n = 1:size(Kusu,1)
%     MC.d13cDIC(n,1:rep) = randn(1,rep) .* Kusu.d13cDIC_std(n) + Kusu.d13cDIC(n);
%     MC.d13cDOC(n,1:rep) = randn(1,rep) .* 0.2 + Kusu.d13cDOC(n);    % the 2 sigma percision of d13C-DOC in uOttawa is 0.4 permil
%     MC.DIC (n,1:rep) = randn(1,rep) .* Kusu.DIC_std(n) + Kusu.DIC(n);
%     MC.TA (n, 1:rep) = randn(1,rep) .* Kusu.TA_std (n) + Kusu.TA(n);
%     MC.DOC (n,1:rep) = randn(1,rep) .* Kusu.DOC_stdev (n) + Kusu.DOC(n);
%     MC.salinity (n,1:rep) = randn(1,rep) .* 0.01 + Kusu.salinity(n);
%     MC.T (n, 1:rep) = randn(1,rep) .* 0 + Kusu.temperature(n);
%     MC.d13c_carb(n,1:rep) = randn(1,rep) .* 0 + 0;
% end
% 
% 
% for n = 1:size(Hantu,1)
%     w = n + size(Kusu,1);
%     MC.d13cDIC(w,1:rep) = randn(1,rep) .* Hantu.d13cDIC_std(n) + Hantu.d13cDIC(n);
%     MC.d13cDOC(w,1:rep) = randn(1,rep) .* 0.2 + Hantu.d13cDOC(n);
%     MC.DIC (w,1:rep) = randn(1,rep) .* Hantu.DIC_std(n) + Hantu.DIC(n);
%     MC.TA (w, 1:rep) = randn(1,rep) .* Hantu.TA_std (n) + Hantu.TA(n);
%     MC.DOC (w,1:rep) = randn(1,rep) .* Hantu.DOC_std (n) + Hantu.DOC(n);
%     MC.salinity (w,1:rep) = randn(1,rep) .* 0.01 + Hantu.salinity(n);
%     MC.T (w, 1:rep) = randn(1,rep) .* 0 + Hantu.temperature(n);
%     MC.d13c_carb(w,1:rep) = randn(1,rep).* 0 + 0;
% end


%% (2) for the actual measurements (use the long-term precision of each measurement)

for n = 1:size(Kusu,1)
    MC.d13cDIC(n,1:rep) = randn(1,rep) .* 0.2 + Kusu.d13cDIC(n);   % the long-term precision in NTU is 0.2 permil
    MC.d13cDOC(n,1:rep) = randn(1,rep) .* 0.2 + Kusu.d13cDOC(n);    % the 2 sigma precision of d13C-DOC in uOttawa is 0.4 permil
    MC.DIC (n,1:rep) = randn(1,rep) .* Kusu.DIC(n) .* 0.0015 + Kusu.DIC(n);  %the precision for DIC CRM is 0.15%
    MC.TA (n, 1:rep) = randn(1,rep) .* Kusu.TA(n) .* 0.0013  + Kusu.TA(n);
    MC.DOC (n,1:rep) = randn(1,rep) .* 3.9 + Kusu.DOC(n);
    MC.salinity (n,1:rep) = randn(1,rep) .* 0.01 + Kusu.salinity(n);
    MC.T (n, 1:rep) = randn(1,rep) .* 0 + Kusu.temperature(n);
    MC.d13c_carb(n,1:rep) = randn(1,rep) .* 0 + 0;
end


for n = 1:size(Hantu,1)
    w = n + size(Kusu,1);
    MC.d13cDIC(w,1:rep) = randn(1,rep) .* 0.2 + Hantu.d13cDIC(n);
    MC.d13cDOC(w,1:rep) = randn(1,rep) .* 0.2 + Hantu.d13cDOC(n);
    MC.DIC (w,1:rep) = randn(1,rep) .* Hantu.DIC(n) .* 0.0015 + Hantu.DIC(n);
    MC.TA (w, 1:rep) = randn(1,rep) .* Hantu.TA(n) .* 0.0013 + Hantu.TA(n);
    MC.DOC (w,1:rep) = randn(1,rep) .* 3.9 + Hantu.DOC(n);
    MC.salinity (w,1:rep) = randn(1,rep) .* 0.01 + Hantu.salinity(n);
    MC.T (w, 1:rep) = randn(1,rep) .* 0 + Hantu.temperature(n);
    MC.d13c_carb(w,1:rep) = randn(1,rep).* 0 + 0;
end

%% calculate the conservative mixing DIC and TA
MC.DIC_cons = DIC_cons_calc(MC.salinity);
MC.TA_cons = TA_cons_calc(MC.salinity);
MC.f_riv = 1 - MC.salinity./MC.sal_mar;
MC.f_mar = 1 - MC.f_riv;
MC.d13c_cons = d13c_mix (MC.dic_riv, MC.d13c_riv, MC.f_riv, MC.dic_mar, MC.d13c_mar, MC.f_mar);

%% calculate the tDOC
MC.R_peat = R_calc(MC.d13c_peat);
MC.R_DOC_mar = R_calc(MC.d13cdoc_mar);   % uncertainty matrix for the R marine endmember of d13c-DOC
MC.R_DOC_meas = R_calc(MC.d13cDOC);  %generate the uncertainty matrix for the meausred R of DOC

MC.tDOC = (MC.R_DOC_meas - MC.R_DOC_mar) .* MC.DOC ./ (MC.R_peat - MC.R_DOC_mar); %calculate the amount of tDOC 
MC.tDOC = unitqlo(MC.tDOC, MC.salinity, MC.T, 1);%convert to umol/kg


%% calculate the carbon budgets, i.e. remineralization and CO2 outgassing
MC.DIC_dev = MC.DIC - MC.DIC_cons;
MC.TA_dev = MC.TA - MC.TA_cons;
MC.e = 23.644 - 9701.5 ./ (MC.T + 273.15); %Rau et al. 1996
MC.a = MC.e ./ 1000 + 1;

A = 1 ./ MC.DIC_cons .* (MC.d13c_tDIC - MC.d13c_cons); 
B = 1 ./ MC.DIC_cons .* 1000 .* (MC.a - 1);
C = 1 ./ MC.DIC_cons .* (MC.d13c_carb - MC.d13c_cons);
        
MC.REM = ((MC.d13cDIC - MC.d13c_cons) - B.*MC.DIC_dev + 0.5 .* MC.TA_dev .* (C - B)) ./ (A - 1.0125 .* B + 0.0125 .* C);
MC.DIS = 0.5 .* MC.TA_dev + 0.0125*MC.REM;
MC.CO2 = MC.REM + MC.DIS - MC.DIC_dev;
MC.REM_remain = MC.REM - MC.CO2;


%% calculate the standard deviations
kk = table();
kk.date = vertcat(Kusu.Date, Hantu.Date);
kk.location = vertcat(Kusu.location, Hantu.location);
kk.tDOC = mean(MC.tDOC, 2);
kk.tDOC_std = std(MC.tDOC, 0, 2);

kk.REM = mean(MC.REM,2, 'omitnan');
kk.REM_std = std(MC.REM, 0 , 2, 'omitnan');

kk.CO2 = mean(MC.CO2, 2, 'omitnan');
kk.CO2_std = std(MC.CO2, 0, 2, 'omitnan');

kk.REM_remain = mean(MC.REM_remain, 2, 'omitnan');
kk.REM_remain_std = std(MC.REM_remain, 0, 2, 'omitnan');

kk.DIS = mean(MC.DIS, 2, 'omitnan');
kk.DIS_std = std(MC.DIS, 0, 2, 'omitnan');

%% calculate the monthly average of the standard deviations
month = datetime('15-Jan-2019');
MC.MthAve = table();

n = 1;
while datenum(month - kk.date(end))<30
    index = find(strcmp(cellstr(datestr(kk.date,'mmm-yyyy')),datestr(month,'mmm-yyyy'))==1);
    MC.MthAve.date(n) = month;
    MC.MthAve{n,2:11} = mean(kk{index,3:end},1,'omitnan');
    
    month = month + (datetime('01-Feb-2019') - datetime('01-Jan-2019'));%loop once, move to the next month
    n = n + 1;
    clear index
end

%change variablenames
for n = 2:size(MC.MthAve,2)
    MC.MthAve.Properties.VariableNames{n} = kk.Properties.VariableNames{n+1};
end




%% compare the result from MC and from the actual measurement
load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\DOC mineralizatoin\tDOC_budget_results.mat','mthave');
mthave = mthave(1:end-1,:);
figure;
plot(mthave.REM,MC.MthAve.REM,'ko');hold on
plot([0 100],[0 100],'k-')

%% arrange the date
for n = 1:size(MC.MthAve,1)
    temp_char = char(MC.MthAve.date(n));
    temp_char = ['15' temp_char(3:end)];
    MC.MthAve.date(n) = cellstr(temp_char);
end


%% draw the bar plot
% the actual monthly mean is calculated from the actual measurements, not
% taken as the results of the MC simulation
figure('color','w');
x = datenum(MC.MthAve.date) - datenum(MC.MthAve.date(1));
y_bar = horzcat(mthave.tDOC, mthave.REM_remain, mthave.CO2);
p = bar(x, y_bar, 0.65, 'stacked', 'facecolor','flat','edgecolor','none');hold on
p(1).FaceColor = [69 137 148]/255;
p(2).FaceColor = [178 200 187]/255;
p(3).FaceColor = [102 179 255]/255;

%only show one-sided error bar for illustration
errorbar(x, mthave.tDOC, zeros(size(x)),MC.MthAve.tDOC_std, 'linestyle', 'none', 'color',[69 137 148]/255); hold on
errorbar(x, mthave.tDOC + mthave.REM_remain, zeros(size(x)),MC.MthAve.REM_remain_std, 'linestyle', 'none','color',[178 200 187]/255);hold on
errorbar(x, mthave.tDOC + mthave.REM_remain + mthave.CO2, zeros(size(x)),MC.MthAve.REM_std, 'linestyle','none', 'color', [102 179 255]/255)

ylabel('terrigenous carbon budget (umol/L)')
%set(gca,'xticklabel',cellstr(datestr(MthAve.date,'mmm')));
legend('tDOC remaining','tDOC remineralized as DIC','tDOC remineralized and degassed');
ylabel('terrestrial carbon budget (umol/L)')
xlabel('Time')
xlim([-30 datenum(datetime('31-Dec-2020') - MC.MthAve.date(1))])
YTICK = [0 20 40 60 80 100 120];
XTICKLABEL = datestr(MC.MthAve.date(1:3:end),'mmm');
XTICK = x(1:3:end);
set(gca,'fontsize',14,'xtick',XTICK,'xticklabel',XTICKLABEL,'ytick',YTICK);


%%
temp = vertcat(Kusu.temperature, Hantu.temperature);
MC.T = repmat(temp, 1, rep);
%% extrapolate the riverine DOC endmember and calculate the uncertainty
% do the regression for each trial
MC.riv_doc = zeros(1,rep);
MC.riv_doc = nan;
for n = 1:rep
    x = MC.salinity(:,n);
    x = horzcat(ones(size(x,1),1),x);
    y = MC.REM(:,n) + MC.tDOC(:,n);
    y = unitqlo(y, MC.salinity(:,n), MC.T(:,n), 2);
    B = regress(y,x);
    MC.riv_doc(n) = B(1);
end

MC.riv_doc_std = std(MC.riv_doc,0,2);


%%
save('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\DOC mineralizatoin\tDOC budget with Monte Carlo.mat');
save('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat','Kusu','Hantu');
%% nested functions
function [d13c_mixed] = d13c_mix (pool_A, d13c_A, f_A, pool_B, d13c_B, f_B)   
%calculate the resulting d13c if two pools of DIC is mixed together
   
    global MC
    R_A = (d13c_A ./ 1000 + 1) .* MC.R_vpdb;
    R_B = (d13c_B ./ 1000 + 1) .* MC.R_vpdb;
    
    C13 = pool_A .* R_A ./(1+R_A) .* f_A + pool_B .* R_B ./(1+R_B) .* f_B;
    C12 = pool_A .* R_A ./(1+R_A) ./ R_A .* f_A + pool_B .* R_B ./(1+R_B) ./ R_B .* f_B;
    
    R_mix = C13 ./ C12;
    d13c_mixed = (R_mix ./ MC.R_vpdb - 1) .* 1000;
end
 
 
function DIC_cons = DIC_cons_calc(salinity)
    global MC
    f_river = 1 - salinity ./ MC.sal_mar;
    f_marine = salinity ./ MC.sal_mar;
    DIC_cons = f_river .* MC.dic_riv + f_marine .* MC.dic_mar;
end
 
function TA_cons = TA_cons_calc (salinity)
    global MC
    f_river = 1 - salinity ./ MC.sal_mar;
    f_marine = salinity ./ MC.sal_mar;
    TA_cons = f_river .* MC.ta_riv + f_marine .* MC.ta_mar;
end


function R = R_calc (d13c)
    global MC
    R = (d13c ./ 1000 + 1) .* MC.R_vpdb;
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


