clear
cd('F:\NTU\Research\Kusu Hantu Biogeochem\SeaFET pH\pH New reprocessed data')
[Num{1},Text{1}]=xlsread('2015-07 to 2015-08_pro.xlsx');
[Num{2},Text{2}]=xlsread('2015-09 to 2015-11_pro.xlsx');
[Num{3},Text{3}]=xlsread('2015-12 to 2016-02_pro.csv');
%[Num{4},Text{4}]=xlsread('2016-03 to 2016-06_pro.csv'); %this file with
%problematic pH value
[Num{4},Text{4}]=xlsread('2016-10 to 2016-12_pro.xlsx',2);  %The first 10 measurements seem to be wrong
[Num{5},Text{5}]=xlsread('2016-12 only_pro.csv'); 
[Num{6},Text{6}]=xlsread('2016-12 to 2017-05_pro.csv');
[Num{7},Text{7}]=xlsread('2017-06 to 2017-08_pro.csv');
[Num{8},Text{8}]=xlsread('2017-08 to 2017-10_pro.csv');
[Num{9},Text{9}]=xlsread('2017-11 to 2017-12_pro.csv');
[Num{10},Text{10}]=xlsread('2017-12 to 2018-01_pro.csv');
[Num{11},Text{11}]=xlsread('2018-01 to 2018-02_pro.csv');
[Num{12},Text{12}]=xlsread('2018-02 to 2018-03_pro.csv');
[Num{13},Text{13}]=xlsread('2018-07 to 2018-08_pro.csv');
[Num{14},Text{14}]=xlsread('2018-08 to 2018-09_pro.csv');

%% 2015-09 to 2015-11 each time point has 5 samples, take the average of them (and also 2015-07 to 2015-09)
%Now the xls file is the correct one
% j=1;
% for i=1:5:40571
%     Num0p(j,:)=mean(Num0(i:i+4,1:9));
%     j=j+1;
% end

%% Uniform the date format and remove the first 20 measurements and last 20 measurements
% File 1-5 have the same format
%pH_Data=struct();

for i=1:3
    Date_numeric = Num{i}(10:end-10,1);
    Time = Num{i}(10:end-10,2);
    %conversion to normal format and local time
    for j=1:length(Date_numeric)
        if Date_numeric(j,1)>2015000 && Date_numeric(j,1)<2016000
            pH_Data.Date{i}(j,1)=datetime('01-Jan-2015')+(Date_numeric(j,1)-2015000)-1;   %convert to dd-mmm-yyyy
            pH_Data.Date{i}(j,1)=pH_Data.Date{i}(j,1)+Time(j,1)*(datetime('01-Jan-2018 01:00:00')-datetime('01-Jan-2018 00:00:00'))+(datetime('01-Jan-2018 08:00:00')-datetime('01-Jan-2018 00:00:00')); %include the hour/min/sec and convert to local time 
           
        elseif Date_numeric(j,1)>2016000
            pH_Data.Date{i}(j,1)=datetime('01-Jan-2016')+(Date_numeric(j,1)-2016000)-1;   %convert to dd-mmm-yyyy
            pH_Data.Date{i}(j,1)=pH_Data.Date{i}(j,1)+Time(j,1)*(datetime('01-Jan-2018 01:00:00')-datetime('01-Jan-2018 00:00:00'))+(datetime('01-Jan-2018 08:00:00')-datetime('01-Jan-2018 00:00:00')); %include the hour/min/sec and convert to local time 
              
        end
    end
    %compile the data
    pH_Data.Internal_pH{i}=Num{i}(10:end-10,6);  %column 6 is internal
    pH_Data.External_pH{i}=Num{i}(10:end-10,7);
    pH_Data.SeaFET_Temperature{i}=Num{i}(10:end-10,5); %record its own temperature
    pH_Data.ID{i}=[datestr(pH_Data.Date{i}(1,1),'mmm-yyyy'),' to ',datestr(pH_Data.Date{i}(end,1),'mmm-yyyy')]
    i
    clear Date_numeric Time 
end

%%
for i=4:5
    Num{i}=Num{i}(11:end,:);
    Date_numeric = Num{i}(10:end-10,1);
    Time = Num{i}(10:end-10,2);
    %conversion to normal format and local time
    for j=1:length(Date_numeric)
        if Date_numeric(j,1)>2015000 && Date_numeric(j,1)<2016000
            pH_Data.Date{i}(j,1)=datetime('01-Jan-2015')+(Date_numeric(j,1)-2015000)-1;   %convert to dd-mmm-yyyy
            pH_Data.Date{i}(j,1)=pH_Data.Date{i}(j,1)+Time(j,1)*(datetime('01-Jan-2018 01:00:00')-datetime('01-Jan-2018 00:00:00'))+(datetime('01-Jan-2018 08:00:00')-datetime('01-Jan-2018 00:00:00')); %include the hour/min/sec and convert to local time 
           
        elseif Date_numeric(j,1)>2016000
            pH_Data.Date{i}(j,1)=datetime('01-Jan-2016')+(Date_numeric(j,1)-2016000)-1;   %convert to dd-mmm-yyyy
            pH_Data.Date{i}(j,1)=pH_Data.Date{i}(j,1)+Time(j,1)*(datetime('01-Jan-2018 01:00:00')-datetime('01-Jan-2018 00:00:00'))+(datetime('01-Jan-2018 08:00:00')-datetime('01-Jan-2018 00:00:00')); %include the hour/min/sec and convert to local time 
              
        end
    end
    %compile the data
    pH_Data.Internal_pH{i}=Num{i}(10:end-10,6);  %column 6 is internal
    pH_Data.External_pH{i}=Num{i}(10:end-10,7);
    pH_Data.SeaFET_Temperature{i}=Num{i}(10:end-10,5); %record its own temperature
    pH_Data.ID{i}=[datestr(pH_Data.Date{i}(1,1),'mmm-yyyy'),' to ',datestr(pH_Data.Date{i}(end,1),'mmm-yyyy')]
    i
    clear Date_numeric Time 
end


%% File No.6 -12 have different time format
for i=6:12
    Num{i}=Num{i}(11:end,:);
    Text{i}=Text{i}(15:end,:);
    pH_Data.Date{i} = datetime(datestr(Text{i}(6:end-5,3)))+(datetime('01-Jan-2018 08:00:00')-datetime('01-Jan-2018 00:00:00'));  %convert to local time    
    pH_Data.Internal_pH{i} = Num{i}(6:end-5,6);
    pH_Data.External_pH{i} = Num{i}(6:end-5,7);
    pH_Data.SeaFET_Temperature{i} = Num{i}(6:end-5,5);   %use its own temperature
    pH_Data.ID{i}=[datestr(pH_Data.Date{i}(1,1),'mmm-yyyy'),' to ',datestr(pH_Data.Date{i}(end,1),'mmm-yyyy')]
    i
end
%%
for i=13:14
    Num{i}=Num{i}(11:end,:);
    Text{i}=Text{i}(15:end,:);
    pH_Data.Date{i} = datetime(datestr(Text{i}(6:end-5,3),0))+(datetime('01-Jan-2018 08:00:00')-datetime('01-Jan-2018 00:00:00'));  %convert to local time    
    pH_Data.Internal_pH{i} = Num{i}(6:end-5,6);
    pH_Data.External_pH{i} = Num{i}(6:end-5,7);
    pH_Data.SeaFET_Temperature{i} = Num{i}(6:end-5,5);   %use its own temperature
    pH_Data.ID{i}=[datestr(pH_Data.Date{i}(1,1),'mmm-yyyy'),' to ',datestr(pH_Data.Date{i}(end,1),'mmm-yyyy')]
    i
end

%% Arrange the dates
pH_Data.Date_all=pH_Data.Date{1};
for i=2:14
    pH_Data.Date_all=vertcat(pH_Data.Date_all,pH_Data.Date{i});
end

%% Compile all the data
pH_Data.Internal_all=pH_Data.Internal_pH{1};
pH_Data.External_all=pH_Data.External_pH{1};
pH_Data.Temperature_all = pH_Data.SeaFET_Temperature{1};
for i=2:14
     pH_Data.Internal_all=vertcat(pH_Data.Internal_all,pH_Data.Internal_pH{i});
     pH_Data.External_all=vertcat(pH_Data.External_all,pH_Data.External_pH{i});
     pH_Data.Temperature_all=vertcat(pH_Data.Temperature_all,pH_Data.SeaFET_Temperature{i});
end

%% Save the data
% xlswrite('pH full record Data(no gaps).xlsx',[Internal_all,External_all]);
% xlswrite('pH full record dates (no gaps).xlsx',cellstr(Dates_all))
cd('F:\NTU\Research\Kusu Hantu Biogeochem\SeaFET pH\pH New reprocessed data\pH Matlab workspace time series')
save SeaFET_Jul-2015_to_Sep_2018.mat

%% Next time starts from here, just add the new data to the full dataset
clear
load ('F:\NTU\Research\Kusu Hantu Biogeochem\SeaFET pH\pH New reprocessed data\pH Matlab workspace time series\SeaFET_201507_to_202007.mat')

%% PLS CREATE A NEW LIST OF UTC DT IN THE CORRECTED PH FILE BEFORE IMPORTING
%import new data
cd('F:\NTU\Research\Kusu Hantu Biogeochem\SeaFET pH\pH New reprocessed data')
[NewNum,NewText]=xlsread('SeaFET Kusu 2020-07 to 2020-10_pro.csv');
if NewNum(12,5)<10
    fprintf('ERROR ERROR ERROR ERROR ERROR ERROR ERROR')
end


%% The imported time is UTC rather than local time, so i have to convert it
% to local time, and then add the new data to the full dataset
newdata = table();
newdata.DT = (datetime(datestr(NewText(16:end-7,3),0))+(datetime('01-Jan-2018 08:00:00')-datetime('01-Jan-2018 00:00:00')));
newdata.int = NewNum(12:end-7,6); %remove the first and last few measurements
newdata.ext = NewNum(12:end-7,7);
newdata.temperature = NewNum(12:end-7,5);

ph = vertcat(ph, newdata);

%% save data
save('F:\NTU\Research\Kusu Hantu Biogeochem\SeaFET pH\pH New reprocessed data\pH Matlab workspace time series\SeaFET_201507_to_2020_10.mat')

%% import the calculated pH
Kusu_raw = readtable ('F:\NTU\Research\Kusu Hantu Biogeochem\Kusu Hantu Biogeochem summary\Kusu Hantu 5m Biogeochem data summary.xlsx','Sheet','Kusu');
Kusu = Kusu_raw(strcmp('Kusu',Kusu_raw.location)==1,:);


%% Plot the full dataset
%% remove abnormal values
ph.int(ph.int>8.1 | ph.int<7.8) = NaN;
ph.ext(ph.ext>8.1 | ph.ext<7.8) = NaN;

%% do the plotting
Xtick = {'15-Jul-2015','15-Jan-2016','15-Jul-2016','15-Jan-2017','15-Jul-2017','15-Jan-2018','15-Jul-2018','15-Jan-2019','15-Jul-2019','15-Oct-2019','15-jul-2020'};
Xtick = datetime(Xtick);

figure('color','w')
h1=plot(datetime('01-Jan-2000'),1,'ro','markerfacecolor','r')
hold on
h2=plot(datetime('01-Jan-2000'),1,'bo','markerfacecolor','b')
hold on
 plot(datetime(dm_ph.DT),dm_ph.int,'r.')
 hold on
plot(datetime(dm_ph.DT),dm_ph.ext,'b.')
hold on
h3=plot(datetime(datestr(Kusu.Date,0)), Kusu.pH_calc,'go','markerfacecolor','g')

xlim([datetime('01-jul-2015') datetime('31-Dec-2020')])
ylim([7.8  8.1])
% set(gca, 'XTick', Xtick)
% set(gca,'Xticklabel',datestr(Xtick,'mmm-yy'))

set(gca,'Fontsize',14)
title('pH record from Jul-2015 to Oct-2020')
ylabel('pH')
xlabel('Date')
legend([h1 h2 h3],{'Internal sensor','External sensor','Calculated pH by CO2SYS'})
%legend([h1 h2],{'Internal sensor','External sensor'})

%% quality control, remove the bad measurements

%% calculate the dailymean pH
new_dm_ph = table();
newdata.DT = cellstr(datestr(newdata.DT,'dd-mmm-yyyy'));

k = 1; i = 1; j = 1;
while k < size(newdata,1)+1
    if k == size(newdata,1) || strcmp(newdata.DT(k),newdata.DT(k+1)) == 0 
        new_dm_ph.DT(j) = newdata.DT(k);
        new_dm_ph.int(j) = mean(newdata.int(i:k),'omitnan');
        new_dm_ph.ext(j) = mean(newdata.ext(i:k),'omitnan');

        j = j + 1;
        i = k + 1;
    end
    k = k + 1;
end

% then add the new dailymean to the old dailymean dataset
load('F:\NTU\Research\Kusu Hantu Biogeochem\SeaFET pH\pH New reprocessed data\dailymean pH\dm_ph.mat','dm_ph')
dm_ph = vertcat(dm_ph,new_dm_ph);
save('F:\NTU\Research\Kusu Hantu Biogeochem\SeaFET pH\pH New reprocessed data\dailymean pH\dm_ph.mat','dm_ph')

%% we can choose the appropriate sensor to be plotted during different periiods
figure('color','w')
plot(datetime(vertcat(pH_Data.Date_all(1:109328),pH_Data.Date_all (121734:end))),vertcat(pH_Data.External_all(1:109328), pH_Data.Internal_all(121734:end)),'b.')
hold on
plot(datetime(datestr(Kusu.Date,0)), Kusu.pH_calc,'ko','markerfacecolor','k')

xlim([datetime('01-jul-2015') datetime('16-Sep-2019')])
%xlim([datetime('01-Oct-2017') datetime('15-Feb-2019')])

ylim([7.75  8.10])
% set(gca, 'XTick', Xtick)
% set(gca,'Xticklabel',datestr(Xtick,'dd-mmm-yy'))

set(gca,'Fontsize',14)
%title('pH record from Jul-2015 to Sep-2019')
ylabel('pH')
xlabel('Date')
legend('pH logger','calculated pH')

%%
figure('color','w')

plot(datetime(pH_Data.Date_all(127777:127920)),pH_Data.External_all(127777:127920),'b.')
