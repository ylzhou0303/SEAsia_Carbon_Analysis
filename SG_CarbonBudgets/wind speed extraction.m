%% EXTRACTION OF WINDSPEED FROM SATELLITE DATA
% In order to calculate the outgassing at Kusu, I extract the windspeed from the CYGNSS satellite data for our fieldwork dates.

%%
cd('F:\NTU\Research\Kusu Hantu Biogeochem\Kusu Hantu Biogeochem summary')
Kusu_raw = readtable ('Kusu Hantu 5m Biogeochem data summary.xlsx','Sheet','Kusu');
Kusu = Kusu_raw(strcmp('Kusu',Kusu_raw.location)==1,:);

% calculate the day of year
m = size(Kusu);
for i = 1:m(1)
    
    Kusu.doy(i) = datenum(Kusu.Date(i)) - datenum(datetime(['01-Jan-' num2str(year(Kusu.Date(i)))])) + 1;
end


%% specify a region, then calculate the mean wind speed in this region
clear joey
cd('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\main figures\F2 pCO2\CYGNSS L2');
[finfo outstrct]=read_nc_file_struct('cyg.ddmi.s20210125-000000-e20210125-235959.l2.wind-mss.a21.d21.nc');

yung = table();
yung.wind_speed = outstrct.wind_speed;
yung.lat = outstrct.lat;
yung.lon = outstrct.lon;
joey = yung (yung.lat > 0.5 & yung.lat < 1.5 & yung.lon> 102.5 & yung.lon < 105.5, :);
meanws = mean(joey.wind_speed, 'omitnan');

% feed the windspeed into Kusu
% Kusu.wind_speed (k) = meanws;
% k = k+1;


%%
clear outstrct finfo

save('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\pCO2\CYGNSS L2\Kusu_windspeed')
save('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat','Kusu','Hantu')
