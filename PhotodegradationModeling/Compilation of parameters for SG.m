% COMPILATION OF ALL THE PARAMETERS NEEDED FOR THE MODELLING: PABS, bb,
% AQY, E(scalar, below surface), DOC, CDOM
% wavelength range: 300 - 700nm (PABS measurement starts from 300nm)
% It means that the light between 290 and 300nm will have to be neglected for the
% modelling

%% import the AQY for CDOM
clear
load('F:\NTU\Research\Photodegradation experiment\photodeg modelling\MeanCDOMAQY.mat','MeanCDOM_AQY','MeanCDOM_AQY_std')

AQY_CDOM_250_700 = table();
AQY_CDOM_250_700.wl = [250:700]';   % model the change in CDOM spectrum of 250-700nm. I can only start from 250nm, because the CDOM data of Sematan, Samunsam starts from 250nm, due to the NaN3.
SolarLight_wl = [300:700];   % For the AQY spectrum, it should match the solar light spectrum
AQY_CDOM_250_700{:,2:402} = MeanCDOM_AQY(21:end,:);


%% import AQY for DOC loss
load('F:\NTU\Research\Photodegradation experiment\photodeg modelling\AQYs_DOC.mat');

%% FOR Singapore water modeling
%% PABS and bb
load('F:\NTU\Research\Photodegradation experiment\Kd\Kdmodelling.mat')
% get the average PABS and bb
meanPABS = table();
meanPABS.wl = [300:700]';
joey = mean(PABS{:,9:409},1);
yung = std(PABS{:,9:409},0,1);
meanPABS.pabs = joey';
meanPABS.stdev = yung';

meanbb = table();
meanbb.wl = [300:700]';
stephanie = mean(bbfit{:,4:404},1);
cheng = std(bbfit{:,4:404},0,1);
meanbb.bb = stephanie';
meanbb.stdev = cheng';

W = W(101:501,:);


%% E (below surface, scalar)
location = 'SG';
condition = 'clear';  % need to change the condition here if needed

matname = ['Spectra_od_',location,'_',condition,'.mat'];
load(['F:\NTU\Research\Photodegradation experiment\TUV model\TUV outputs\',matname],'Spectra_od','SZA');

SolarSpc = table();
SolarSpc.UTC = SZA.UTC;
SolarSpc.SZA = SZA.sza;

Spectra_tailored = Spectra_od(:,71:471);

% convert to Einst
plank = 6.62607004*10^-34; %Plank constant
c = 299792458;  % speed of light   m/s
avog = 6.022*10^23; % Avogadro constant
wl = [300:700];

Spc_NumPho = Spectra_tailored .* ( wl .* 10^-9) ./(plank * c);  %convert energy to number of photons   n m-2 s-1
Spc_Einst = Spc_NumPho ./ avog;       % convert from number of photons to mol of photons   mol m-2 s-1
SolarSpc{:,3:403} = Spc_Einst;  % scalar solar spectrum below surface, unit: mol photon m-2 nm-1 s-1, range: 300-700 nm

%% FOR SG, CALCULATE THE DOC AND CDOM AT T0, assuming conservative mixing at salinity of 29
DOC_river = 890;  %unit: umol L-1
DOC_t0 = DOC_river .* (1 - 29/33);

load('F:\NTU\Research\Photodegradation experiment\Experimental Data\Maludam photodegradation May-2018\AQYs.mat','MalPhbl');
CDOM_t0 = table();
CDOM_t0.wl = [250:700]';
CDOM_t0.river = MalPhbl.a_raw{1,33:483}' ./ 3250 .* 890;
CDOM_t0.shelf = CDOM_t0.river .* (1-29/33) ;  % normalize the CDOM to the level of Sumatra average


%% water depth
MLD = 19.66;
MLD_std = 12.79;

%%
savename = ['PARs_',location,'_',condition,'.mat'];
clearvars -except AQY_CDOM_250_700 AQYs_DOC SolarSpc DOC_t0 CDOM_t0 meanPABS meanbb W MLD MLD_std DOC_river MeanCDOM_AQY MeanCDOM_AQY_std savename
save(['F:\NTU\Research\Photodegradation experiment\photodeg modelling\',savename]);
