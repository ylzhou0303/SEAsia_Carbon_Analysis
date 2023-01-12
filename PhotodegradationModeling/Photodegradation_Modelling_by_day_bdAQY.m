% CALCULATION BY DAY
% GENERATE A MATRIX OF MEAN ABSORBED PHOTON FOR EACH DAY, ONLY CALCULATE
% THE CHANGES IN DOC AND CDOM BY DAY, NOT BY HOUR

% photodegradation modelling

%% import data
clear
location = 'Talang';  %SG or Talang
condition = 'cloud';
season = 'winter';  % for Talang, enter spring, summer, fall or winter; for SG, just leave it empty

matname = ['PARs_',location,'_',condition,'.mat'];
load (['C:\NTU\Research\Photodegradation experiment\photodeg modelling\',matname]);

% SolarSpc: scalar downward irrandiance just below the water surface, unit:
% mol m-2 nm-1 s-1

%% set up parameters
if strcmp(location,'Talang')
   days = 90;  %modelling period
elseif strcmp(location,'SG')
    days = 730;
end

S = 1;  %surface area, 1 m2
MLD = MLD;  % water depth
V = S * MLD; % volume of the water cubic we want to model
SolarSpc_wl = 300:700; %sunlight spectrum wavelength
CDOM_wl = 250:700; %CDOM wavelength

E_all = reshape(SolarSpc{:,3:end}', 401, 24, 12);  
%one page corresponds to one month, in each page, each COLUMN corresponds to each hour, each ROW corresponds to each wavelength
SZA_all = reshape(SolarSpc.SZA, 1,24,12);

% set up the AQY for DOC and CDOM    
AQY_DOC_spc = AQYs_DOC.Mean';
AQY_CDOM_spc = AQY_CDOM_250_700{:,2:end};

% set up start point for DOC and the CDOM spectrum
DOC_t = table();
DOC_t.DOC_mol = ones(days,1)*nan;
DOC_t.DOC_mol(1) = (DOC_t0 .* 10^-6 .* 10^3) .* V;  %converted from umol/L to mol/m3, and then multiplied by the volume of the water column
DOC_t.DOCloss = ones(days,1)*nan; %to record the DOC loss for each day

CDOM_t = ones(days, 451).*nan;
CDOM_t(1,:) = CDOM_t0.shelf';



%water optical variables
CDOMforKd = CDOM_t0.shelf(51:end);  % unit: m-1 %only use the 300-700nm of CDOM for Kd calculation
%The Kd for underwater light field does not change over time, because it
%has other sources and fresh input of CDOM, but in the modeling, the CDOM
%of the original tDOC fraction would decrease over time due to
%photobleaching
BB = meanbb.bb' + W.bb'; % set up the backscaterring, which does not change over time


% a350_loss = ones(days,1)*nan;
% f_cdom = ones(days,1)*nan;
% doc_loss = ones(days,1)*nan;
% qperday = ones(days,401)*nan;


%% broadband AQY
% use the mean of the broadband AQY calculated from Exp 1-4
AQY_DOC_bd = mean([68.61 69.27 33.34 73] * 10^-6);  %unit: mol C loss/ mol photons

%%
startdate = joey(season); %work out the start date based on the season first (only for Talang region)

for d = 2:days  
    date = (d-1) + startdate;
    month = str2num(datestr(date,'mm'));

    % pick the appropriate solar spectra according to which month you are
    % in
    E = E_all(:,:,month)';
    SZA = SZA_all(:,:,month)';
    
    % calculate the Kd
    CDOMforKd = CDOM_t(d-1, 51:end);  % update the Kd using the CDOM at each day
    A = CDOMforKd + meanPABS.pabs' + W.abs';  % the total absorption is the sum of CDOM, PABS and water abs
    Kd = (1 + 0.005 .* SZA) .* A + 4.18 .* (1 - 0.52 .* exp(-10.8 * A)) .* BB; %Lee et al.2005

    %calculate the fraction of CDOM in the total Kd
    F_CDOM = ((1 + 0.005 .* SZA) .* CDOM_t(d-1,51:end)) ./ Kd; %only take 300-700nm of the CDOM_t to match the length of Kd
    
    % generate the matrix of absorbed photons that are absorbed by CDOM
    Q = E .* S .* (1 - exp(-Kd .* MLD)) .* F_CDOM; %unit: mol nm-1 s-1, absorbed photons per second for each hour
    Qperday = sum(Q .* 3600);  %total amount of photon absorbed per day, unit: mol nm-1 day-1
    
    % multiply the broadband AQY by the absorbed photon and calculate the loss of
    % DOC
    DOC_t.DOCloss(d-1) = sum(Qperday .* AQY_DOC_bd); %unit: mol
    DOC_t.DOC_mol(d) = DOC_t.DOC_mol(d-1) - DOC_t.DOCloss(d-1); %DOC after photodeg at the end of the day, i.e. start of next day
    
    % calculate the CDOM
    cdomloss = sum(Qperday .* AQY_CDOM_spc, 2); %unit: (L m-1)
    cdomloss = cdomloss' ./ (V .* 10^3); %convert the unit from (L m-1) to m-1
    CDOMloss(d-1,:) = cdomloss;
    CDOM_t(d,:) = CDOM_t(d-1,:) - CDOMloss(d-1,:); %calculate the CDOM at next time point by substracting the CDOMloss from the CDOM at the last time point
    
%     a350_loss(d) = CDOMloss_pervol(101);  %350nm
     f_cdom(d-1) = sum(Qperday)./(sum(sum(E.*3600))); %300 - 700nm
%     E_day = sum(E_month .* 3600);
%     f_a350(d) = Qperday(51) ./ E_day(51);
     doc_loss(d-1) = DOC_t.DOCloss(d-1) .* 10^6;  %umol m-2 day-1
%     qperday(d,1:401) = Qperday;
    
end

DOC_t.DOC_umolL = (DOC_t.DOC_mol .* 10^6) ./ (V .* 10^3);%convert from mol to umol/L

%%
savename = ['PhotoDegModel_',location,'_',condition,'_',season,'.mat'];
save(['C:\NTU\Research\Photodegradation experiment\photodeg modelling\model results\BroadbandAQY\',savename]);

%%
figure;
plot([250:700],DOC_t{1:50:end,7:end});hold on
plot([250:700],CDOM_t0.shelf,'k-')

%% nested functions
function [startdate] = joey(season)
    if strcmp(season,'spring')
        startdate = datetime('01-Feb-2019')
    elseif strcmp(season,'summer')
        startdate = datetime('01-May-2019')
    elseif strcmp(season,'fall')
        startdate = datetime('01-Aug-2019');
    elseif strcmp(season,'winter')
        startdate = datetime('01-Nov-2019')
    elseif isempty(season)
        startdate = datetime('01-Jan-2019')
    end
end