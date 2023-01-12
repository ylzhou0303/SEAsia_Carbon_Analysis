cd('F:\NTU\Research\Kusu Hantu Biogeochem\Seabird Sal\processed data for pH correction')
[Num{1},Text{1}]=xlsread('Seabird CTD data_all.xlsx',1);   %create a cell array
[Num{2},Text{2}]=xlsread('Seabird CTD data_all.xlsx',2);
[Num{3},Text{3}]=xlsread('Seabird CTD data_all.xlsx',3);
[Num{4},Text{4}]=xlsread('Seabird CTD data_all.xlsx',4);
[Num{5},Text{5}]=xlsread('Seabird CTD data_all.xlsx',5);  
[Num{6},Text{6}]=xlsread('Seabird CTD data_all.xlsx',6); 
[Num{7},Text{7}]=xlsread('Seabird CTD data_all.xlsx',7);
[Num{8},Text{8}]=xlsread('Seabird CTD data_all.xlsx',8);
[Num{9},Text{9}]=xlsread('Seabird CTD data_all.xlsx',9);
[Num{10},Text{10}]=xlsread('Seabird CTD data_all.xlsx',10);      %deployed on Dec 15, stopped at Dec-23
[Num{11},Text{11}]=xlsread('Seabird CTD data_all.xlsx',11);    %deployed on Jan-12, stopped at Jan-16
[Num{12},Text{12}]=xlsread('Seabird CTD data_all.xlsx',12);    %deployed on Feb-24, stopped on Feb-28
[Num{13},Text{13}]=xlsread('Seabird CTD data_all.xlsx',13);

%% uniform the format of all the dates and remove the first and the last 20 data points
for i=1:13
    Dates{i}=cellstr(datestr(Text{i}(20:end-20),0));
    Salinity{i}=Num{i}(20:end-20,2);
    Temperature{i}=Num{i}(20:end-20,1);
end


%% 2015-09 to 2015-11 each time point has 5 samples, take the average of them (and also 2015-07 to 2015-09)
% j=1;
% for i=1:5:40571
%     Num0p(j,:)=mean(Num0(i:i+4,1:9));
%     j=j+1;
% end

%% Arrange the dates
Dates_all=Dates{1};
for i=2:13
    Dates_all=vertcat(Dates_all,Dates{i});
end

%% Compile all the data
Salinity_all=Salinity{1};
Temperature_all = Temperature{1};
for i=2:13
    Cleaningdate(i-1)=Dates{i-1}(end);
    Salinity_all=vertcat(Salinity_all,Salinity{i});
    Temperature_all=vertcat(Temperature_all,Temperature{i});
end

cd('F:\NTU\Research\Kusu Hantu Biogeochem\Seabird Sal\processed data for pH correction')
save salinity_record_Jul-2015_to_May-2018; 

%% Next time can start from here
% import previous data
clear
load('F:\NTU\Research\Kusu Hantu Biogeochem\Seabird Sal\processed data for pH correction\Matlab workspace time series\Seabird CTD 201507 to 202007.mat')

%% import new data collected from the last field trip
newdata=readtable('F:\NTU\Research\Kusu Hantu Biogeochem\Seabird Sal\processed data for pH correction\Seabird CTD data_compiled_by_month.xlsx','sheet','2020-07 to 2020-10');
newdata.DT = newdata.DT + (datetime('01-jan-2020 08:00:00') - datetime('01-jan-2020 00:00:00'));
saldata = vertcat(saldata,newdata);

%% save data
save('F:\NTU\Research\Kusu Hantu Biogeochem\Seabird Sal\processed data for pH correction\Matlab workspace time series\Seabird CTD 2015_07 to 2020_10.mat')

%% import Valeport CTD data
Kusu_raw = readtable ('F:\NTU\Research\Kusu Hantu Biogeochem\Kusu Hantu Biogeochem summary\Kusu Hantu 5m Biogeochem data summary.xlsx','Sheet','Kusu');
Kusu = Kusu_raw(strcmp('Kusu',Kusu_raw.location)==1,:);


%% Plot the logger data with the CTD data
Xtick = {'15-Jul-2015','15-Jan-2016','15-Jul-2016','15-Jan-2017','15-Jul-2017','15-Jan-2018','15-Jul-2018','15-Jan-2019','15-Jul-2019','15-Oct-2019','15-jul-2020'};
Xtick = datetime(Xtick);
figure('color','w')

h1=plot(datetime('1-Jul-2000'),30,'bo','markerfacecolor','b')
hold on


% highlight the monsoon seasons
% hold on
% for year = 2015:2019 
%     sw_start = datetime(horzcat('01-Jun-',num2str(year)));
%     sw_end = datetime(horzcat('30-Sep-',num2str(year)));
%     ne_start = datetime(horzcat('15-Nov-',num2str(year)));
%     ne_start = datetime(horzcat('28-Feb-',num2str(year+1)));
%     
%     rectangle('position',[datenum(sw_start - datetime('1-Jul-2000')),28,datenum(sw_end - sw_start),34-28],'facecolor',[0.9 0.9 0.9]);
%     hold on
% end
% 
% 

plot(datetime(saldata.DT),saldata.salinity,'.')
hold on
h2=plot(datetime(datestr(Kusu.Date,0)),Kusu.S_Valeport,'go','MarkerFaceColor','g','markersize',8) %Data from Valeport CTD
% set(gca, 'XTick',Xtick)
% set(gca,'Xticklabel',datestr(Xtick,'mmm-yy'))
xlim([datetime('1-Jul-2015') datetime('31-Dec-2020')])

set(gca,'Fontsize',18)
title('Salinity Record at Kusu (Jul-2015 to Oct-2020)')



ylim([28 34])
ylabel('Salinity')
xlabel('Date')
legend([h1 h2],{'Seabird CTD','Valeport CTD'})


%% if the data look good, calculate the dailymean for the newdata
new_dm_sal = table();
newdata.DT = cellstr(datestr(newdata.DT,'dd-mmm-yyyy'));

k = 1; i = 1; j = 1;
while k < size(newdata,1)+1
    if k == size(newdata,1) || strcmp(newdata.DT(k),newdata.DT(k+1)) == 0 
        new_dm_sal.DT(j) = newdata.DT(k);
        new_dm_sal.temperature(j) = mean(newdata.temperature(i:k),'omitnan');
        new_dm_sal.salinity(j) = mean(newdata.salinity(i:k),'omitnan');
        
        j = j + 1;
        i = k + 1;
    end
    k = k + 1;
end
    
%% import old data and add the new data
load ('F:\NTU\Research\Kusu Hantu Biogeochem\Seabird Sal\processed data for pH correction\daily mean salinity\dailymean_sal.mat')
dm_sal = vertcat(dm_sal,new_dm_sal);
save ('F:\NTU\Research\Kusu Hantu Biogeochem\Seabird Sal\processed data for pH correction\daily mean salinity\dailymean_sal.mat','dm_sal')


%% plot with pCO2 records
 Dates_all2 = datetime(Dates_all)
Xtick2=Xtick;

cd('F:\NTU\Research\Kusu Hantu Biogeochem\SAMI pCO2 sensor\processed data');
 load pCO2_records_Mar-2018_to_June-2018.mat

 %%
[AX,h1,h2]=plotyy(Dates_all2(91573:end),Salinity_all(91573:end),Dates_all,Data_all);
set(h1,'marker','.')
set(h1,'linestyle','none')
set(get(AX(1),'ylabel'),'string','Salinity')
set(get(AX(2),'YLABEL'),'STRING','pCO_2 (uatm)')
xlabel('Date')
set(gca,'Xtick',Xtick)
set(gca,'Xticklabel',datestr(Xtick,'dd-mmm-yyyy'))
title('Salinity and pCO_2 records at Kusu')
legend([h1 h2],'Salinity','pCO_2','location','northwest')
set(gca,'fontsize',14)
set(get(AX(2),'YLABEL'),'fontsize',20)

%% find out the salinity of the dates when we do not have the Valeport CTD data
j=1;
for i = 120550:139092
    if strcmp(Text_all(i),'21-Mar-2018')==1
        A2(j)=Data_all(i);
        %B(j)=i;
        j = j+1;
    end
end
if j==1
    fprintf('\nNot found')
end
%salinity=mean(A);   