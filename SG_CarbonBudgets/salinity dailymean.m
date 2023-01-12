
%% SALINITY DAILY MEAN
%import new salinity data to calculate the daily mean
clear
data = readtable('F:\NTU\Research\Kusu Hantu Biogeochem\Seabird Sal\processed data for pH correction\Seabird CTD data_all (use local time here).xlsx','sheet','2019-12 to 2020-02');

%% dailymean of salinity
k=1;
j=1;
data.DT(end+1,1)=cellstr('01-Jan-1900');
for i=1:length(data)
    temp_salinity(j)=data.salinity(i);
    j=j+1;
    if strcmp(datestr(data.DT(i,1),'dd-mmm-yyyy'),datestr(data.DT(i+1),'dd-mmm-yyyy'))==0
        dailymean_sal(k,1)=mean(temp_salinity,'omitnan');
        dailymeandate(k,1)=cellstr(datestr(data.DT(i,1),'dd-mmm-yyyy'));
        k=k+1;
        j=1;
        temp_salinity=[];
    end
end


%% dailymean of temperature
k=1;
j=1;
Date(end+1,1)=cellstr('01-Jan-1900');
for i=1:length(data)
    temp_T(j)=data.temperature(i);
    j=j+1;
    if strcmp(datestr(data.DT(i,1),'dd-mmm-yyyy'),datestr(data.DT(i+1),'dd-mmm-yyyy'))==0
        dailymean_T(k,1)=mean(temp_T,'omitnan');
        dailymeandate(k,1)=cellstr(datestr(Date(i,1),'dd-mmm-yyyy'));
        k=k+1;
        j=1;
        temp_T=[];
    end
end

dm = table();
dm.DT = dailymeandate;
dm.salinity = dailymean_sal;
dm.temperature = dailymean_T;

%% load the daily mean salinity
clear
cd('F:\NTU\Research\Kusu Hantu Biogeochem\Seabird Sal\processed data for pH correction\daily mean salinity')
dmsal = readtable('daily mean salinity.xlsx');
dmsal.salinity(807:837) = NaN;  % these data are believed to be wrong due to large discrepancy with the Valeport sensor

%% add the valeport data
cd('F:\NTU\Research\Kusu Hantu Biogeochem\Kusu Hantu Biogeochem summary')
Kusu_raw = readtable ('Kusu Hantu 5m Biogeochem data summary.xlsx','Sheet','Kusu');
Kusu = Kusu_raw(strcmp('Kusu',Kusu_raw.location)==1,:);

%% plot
% Xtick = {'15-Jul-2015','15-Jan-2016','15-Jul-2016','15-Jan-2017','15-Jul-2017','15-Jan-2018','15-Jul-2018','15-Jan-2019','15-Jul-2019','15-jan-2020'};
% Xtick = datetime(Xtick);

figure('color','w')
plot(datetime(dm_sal.DT(:,1)),dm_sal.salinity(:,1),'bo','markerfacecolor','b','markersize',8);
hold on
plot(datetime(datestr(Kusu.Date,0)),Kusu.S_Valeport,'go','markerfacecolor','g','markersize',8);

% title('Salinity daily mean record of Kusu, Singapore')
xlabel('Time')
ylabel('Salinity')
xlim ([datetime('01-Jul-2015') datetime('31-Dec-2020')])
ylim ()
% legend('Seabird CTD','Valeport CTD')
%set(gca,'xtick',datetime(Xtick),'xticklabel',datestr(Xtick,'mmm-yyyy'),'fontsize',18);


%% feed the daily mean salinity data to the summary file
clear
cd('F:\NTU\Research\Kusu Hantu Biogeochem\SeaFET pH\pH New reprocessed data\dailymean pH')
[Num,Text]=xlsread('dailymean_pH.xlsx');

%%
cd('F:\NTU\Research\Kusu Hantu Biogeochem')
[Num_2,Text_2]=xlsread('Kusu 5m Biogeochem data summary.xlsx');
finddate = cellstr(Text_2(2:end,1));
%%
Date = Text(2:end,1);

for j=1:length(finddate)
    for i=1:length(Date)
        if strcmp(datestr(finddate(j),'dd-mmm-yyyy'),datestr(Date(i,1),'dd-mmm-yyyy'))==1
            temp(j,1:2)=Num(i,1:2);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pCO2
clear
cd('F:\NTU\Research\Kusu Hantu Biogeochem\SAMI pCO2 sensor\processed data\pCO2 Matlab workspace time series')
load pCO2_records_Mar-2018_to_Mar-2019

%%
k=1;
j=1;
Dates_all(end+1,1)=cellstr('01-Jan-1900');
for i=1:length(Data_all)
    temp_pCO2(j)=Data_all(i,1);
    j=j+1;
    if strcmp(datestr(Dates_all(i,1),'dd-mmm-yyyy'),datestr(Dates_all(i+1),'dd-mmm-yyyy'))==0
        dailymeanpCO2(k,1)=mean(temp_pCO2,'omitnan');
        dailymeandatepCO2(k,1)=cellstr(datestr(Dates_all(i,1),'dd-mmm-yyyy'));
        k=k+1;
        j=1;
        temp_pCO2=[];
    end
end

%%
dailymean = dailymeanpCO2;
date = dailymeandatepCO2;
pCO2_dailymean = table(date,dailymean);
save dailymeanpCO2 pCO2_dailymean


%%
load dailymeanpCO2 pCO2_dailymean
%% import the calculated pCO2
cd('F:\NTU\Research\Kusu Hantu Biogeochem\Kusu Hantu Biogeochem summary')
load Kusu

%% remove the values that are too high
f = @pCO2_ErrorDetector
dailymeanpCO2 = arrayfun (f, dailymeanpCO2);

%% plot
Xtick = {'01-Jul-2015','01-Oct-2015','01-Jan-2016','01-Apr-2016','01-Jul-2016','01-Oct-2016','01-Jan-2017','01-Apr-2017','01-Jul-2017','01-Oct-2017','01-Jan-2018','01-Apr-2018','01-Jul-2018','01-Oct-2018','01-Feb-2019'};
Xtick = datetime(Xtick);

figure('color','w')
plot(datetime(dailymeandatepCO2(:,1)),dailymeanpCO2(:,1),'bo','markerfacecolor','b','markersize',10);
hold on
plot(Kusu.Date,Kusu.cal_pCO2,'ro','markerfacecolor','r','markersize',10);

title('pCO2 record from Kusu Island, Singapore (dailymean)')
xlabel('Time')
ylabel('pCO2 (uatm)')
xlim ([datetime('01-Oct-2017') datetime('15-Feb-2019')])
legend ('SAMI pCO2 logger','calculated pCO2')
set(gca,'xtick',datetime(Xtick),'xticklabel',datestr(Xtick,'mmm-yyyy'),'fontsize',14);


%% nest function


function [output] = pCO2_ErrorDetector (x)
    if x>600
        output = NaN;
    else
        output = x;
    end
end
