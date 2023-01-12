% In this file, I use the irradiance measured by radiometer, the AQY of
% Salicylic acid, and Nitrite absorption to predict the production of
% salicylic acid, and compare with the measured data


% import data
% import AQY and absorption
AQY = readtable('C:\NTU\Research\Photodegradation experiment\manuscript\Frontiers in Marine Science\Review and Revisions\Actinometry.xlsx','sheet','Processed AQY');
Abs = readtable('C:\NTU\Research\Photodegradation experiment\manuscript\Frontiers in Marine Science\Review and Revisions\Actinometry.xlsx','sheet','Processed absorption');

%import solar spectrum of the solar simulator
load ('C:\NTU\Research\Photodegradation experiment\Solar simulator irradiance spectrum\solar simulator spectrum.mat','SuntestSpc') % the light spectrum measured by the "old" TriOS radiometer
load('C:\NTU\Research\Photodegradation experiment\Solar simulator irradiance spectrum\new lab radiometer\Suntest TS data\NewRadiometer_SuntestSpc.mat','NewSpcInt');  % the light spectrum measured by the "new" FLAME radiometer

%import water absorption spectrum
load('C:\NTU\Research\Photodegradation experiment\Kd\Kdmodelling.mat','W'); %water absorption

%% calculate the amount of photons that's absorbed by the sample (throughout the entire sample cell) 
V = 30 .* 10^-6;  % volume of sample water, m3
S = (47/2*0.001)^2 * 3.1415926; %area of the cell, unit: m2
PL = V/S; %pathlength of the sample water, m

Conc_nitrite = 0.333 ./ 10^3; % nitrite concentration, 0.333 mM, converted to M
CDOM_Napierian_a = Abs.Abs(1:24)' .* 2.303 .* 100 .* Conc_nitrite;   % the unit of Abs was (mol/L)-1 cm-1, here need to convert to m-1
Kd = W.abs(96:5:211)' + W.bb(96:5:211)' + CDOM_Napierian_a;  % the Kd is the sum of water absorption, backscattering and CDOM
F_CDOM = (1 - exp(-Kd .* PL)) .*  (CDOM_Napierian_a ./ Kd);   % the fraction of irradiance absorbed by CDOM across the entire water column


%% calculate how much photons were absorbed
Q = SuntestSpc.Einst(6:5:121)' .* 0.95 .* S .* F_CDOM;      % calculate using the "old" spectrum
%Q = NewSpcInt.Einst(6:5:121)' .* 0.95 .* S .* F_CDOM;      % calculate
%using the "new spectrum"

%unit: mol photons s-1
% the transmission of light through the quartz material is accounted for by
% a multiplying a factor of 0.95. The material of our quartz cells is
% Spectrosil® Quartz. The transmission of light is 95%. (https://www.starna.com/cells/cell-specifications)


%% now use the AQY to predict the production of salicylic acid
time = 10 * 60; %run exp for 10 minutes
Yield = trapz(AQY.Wavelength(1:24)', Q .* AQY.AQY(1:24)') .* time;  %unit: mol, % 

SA = (Yield .* 10^9) ./ (V .* 10^3); %converted to unit: nM
%%
SA_320 = SA;
%% 
figure;
plot(SuntestSpc.wl, SuntestSpc.Einst)


%% use the time series spectra measured by the new radiometer to calculate the SA production
load('C:\NTU\Research\Photodegradation experiment\Solar simulator irradiance spectrum\new lab radiometer\Suntest TS data\NewRadiometer_SuntestSpc.mat', 'AllSpc_Einst')
E = AllSpc_Einst(6:5:121,:)';
Q = E .* 0.95 .* S .* F_CDOM;  


time = 10; % time interval between each measurement of light spectrum is 10s
Yield = sum(sum(Q .* 5 .* AQY.AQY(1:24)')) .* time;  %unit: mol, % 

SA = (Yield .* 10^9) ./ (V .* 10^3); %converted to unit: nM


%%
save('C:\NTU\Research\Photodegradation experiment\manuscript\Frontiers in Marine Science\Review and Revisions\Actinometry.mat')



