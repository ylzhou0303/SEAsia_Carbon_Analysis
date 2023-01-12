% pCO2 dailymean
% get old dataset
clear
load ('F:\NTU\Research\Kusu Hantu Biogeochem\SAMI pCO2 sensor\processed data\pCO2 dailymean\dailymean_co2', 'dailymean_co2')

% import new data
cd('F:\NTU\Research\Kusu Hantu Biogeochem\SAMI pCO2 sensor\processed data');
data_co2 = readtable ('SAMI data collection_parsed.xlsx', 'sheet','2019-05 to 2019-06');

data_co2.date = datestr(data_co2.localDT, 'dd-mmm-yyyy');
data_co2.date = datetime(data_co2.date);
joey = data_co2.date(1);
dailymean_co2_new = table();

%%
clear temp_date temp_co2 temp_T
k = 1;
while ~(joey > data_co2.date(end))
    temp_date(k) = joey;
    temp_co2 (k) = mean(data_co2.CO2(strcmp(cellstr(data_co2.date), datestr(joey, 'dd-mmm-yyyy')) == 1), 'omitnan');
    temp_T (k) = mean(data_co2.TemperatureC(strcmp(cellstr(data_co2.date), datestr(joey, 'dd-mmm-yyyy')) == 1), 'omitnan');
    joey = joey + (datetime('02-jan-2019') - datetime('01-jan-2019'));
    k = k + 1;    
end

dailymean_co2_new.date = temp_date';
dailymean_co2_new.co2 = temp_co2';
dailymean_co2_new.T = temp_T';

dailymean_co2 = vertcat(dailymean_co2,dailymean_co2_new);


%% plot the daily mean and calculated pCO2
figure('color','w');
plot(dailymean_co2.date, dailymean_co2.co2, 'o', 'markeredgecolor', [102 179 255]/255,'markerfacecolor', [102 179 255]/255, 'markersize', 8);hold on
plot(Kusu.Date, Kusu.pCO2_calc, 'o', 'markeredgecolor',[255 124 124]/255, 'markerfacecolor', [255 124 124]/255, 'markersize', 9)


legend('measured pCO2', 'calculated pCO2')
xlabel('Time')
ylabel('pCO2')
xlim([datetime('01-Oct-2017') datetime('31-Dec-2020')])
ylim([350 800])
XTICK = datetime({'15-oct-2017','15-jan-2018','15-apr-2018','15-jul-2018','15-oct-2018','15-jan-2019','15-apr-2019','15-jul-2019','15-oct-2019','15-jan-2020','15-apr-2020','15-jul-2020','15-oct-2020'});
YTICK = [400 500 600 700 800];
set(gca,'xtick',XTICK, 'xticklabel',datestr(XTICK,'mmm-yy'),'ytick',YTICK)
set(gca,'fontsize',14)

%% save data
cd('F:\NTU\Research\Kusu Hantu Biogeochem\SAMI pCO2 sensor\processed data\pCO2 dailymean')
load dailymean_co2.mat dailymean_co2