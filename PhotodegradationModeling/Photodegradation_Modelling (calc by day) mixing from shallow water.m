% CALCULATION BY DAY
% GENERATE A MATRIX OF MEAN ABSORBED PHOTON FOR EACH DAY, ONLY CALCULATE
% THE CHANGES IN DOC AND CDOM BY DAY, NOT BY HOUR

% photodegradation modelling
load ('F:\NTU\Research\Photodegradation experiment\photodeg modelling\PARs_SG_clear.mat');

% SolarSpc: scalar downward irrandiance just below the water surface, unit:
% mol m-2 nm-1 s-1


%% model the mixing from shallow water to deep water first
days = 720;  %modelling period
S = 1;
depth = 3;
V = S * depth;
SolarSpc_wl = [300:700];
CDOM_wl = [250:700];
salinity = 1;
sal_mar = 33;
DOC = DOC_river .* (1 - salinity / sal_mar); %the DOC at t0
DOC = DOC ./ 10^6 .* 10^3; %convert to mol/m3
DOC_amount = DOC .* S .* depth;

CDOM = CDOM_t0.river' .* (1 - salinity / sal_mar);  % unit: m-1

% set up the solar spectra
E = reshape(SolarSpc{:,3:end}', 401, 24, 12);  
%one page corresponds to one month, in each page, each COLUMN corresponds to each hour, each ROW corresponds to each wavelength
SZA = reshape(SolarSpc.SZA, 1,24,12);

% set up the AQY    
AQY_DOC_spc = AQY_DOC.mean(3:end)';
AQY_CDOM_spc = AQY_CDOM_250_700{:,2:end};


DOC_change = table();
DOC_change.time = [1:days]';
DOC_change.depth = ones(days,1)*nan;
DOC_change.salinity = ones(days,1)*nan;
DOC_change.DOC = ones(days,1)*nan;
DOC_change.DOC_umolL = ones(days,1)*nan;
DOC_change.DOC_amount = ones(days,1)*nan;


%%
d = 1;
while ~ (depth > 24)
    date = (d-1) + datetime('01-Jan-2019');
    month = str2num(datestr(date,'mm'));

    % pick the appropriate solar spectra
    E_month = E(:,:,month)';   % the hourly solar spectra for this month
    SZA_month = SZA(:,:,month)';  % the hourly SZA for this month
    
    % calculate the Kd
    CDOMforKd = CDOM(1,51:end); %only use the 300-700nm of CDOM for Kd calculation
    A = CDOMforKd + meanPABS.pabs' + W.abs';
    BB = meanbb.bb' + W.bb';
    Kd = (1 + 0.005 .* SZA_month) .* A + 4.18 .* (1 - 0.52 .* exp(-10.8 * A)) .* BB; %Lee et al.2005
    
    F_CDOM = CDOMforKd ./ Kd;
    % generate the matrix of absorbed photon
    Q = E_month .* S .* (1 - exp(-Kd .* depth)) .* F_CDOM; %unit: mol nm-1 s-1
    Qperday = sum(Q .* 3600);  %total amount of photon absorbed per day, unit: mol nm-1 day-1
    
    % multiply the AQY by the absorbed photon and calculate the loss of
    % DOC
    DOC_amount_loss = sum(Qperday .* AQY_DOC_spc); %mol
    DOC_amount = DOC_amount - DOC_amount_loss;
    DOC = DOC - DOC_amount_loss ./ V; % mol/m3
    DOC_change.DOC(d) = DOC;
    DOC_change.DOC_umolL(d) = DOC .* 10^6 ./ 10^3;
    DOC_change.DOC_amount(d) = DOC_amount;
    
    % the loss of CDOM
    CDOMloss = sum(Qperday .* AQY_CDOM_spc,2); % m-1 L
    CDOMloss_pervol = CDOMloss ./ (V .* 10^3); %unit: m-1
    CDOM = CDOM - CDOMloss_pervol';
    DOC_change{d,7:457} = CDOM;
    
    % salinity
    DOC_change.salinity(d) = salinity;
    
    % now calculate the variables for day 2
    depth_old = depth;
    depth_next = depth_old + 1;  %add one meter to the depth each day, until the depth becomes 24m
    depth = depth_next;
    DOC_change.depth(d) = depth_old;
    
    V_next = S .* depth_next;
    V = V_next;
    
    DOC_next = DOC_amount ./ (V_next);  %mol/m3
    DOC = DOC_next;
    CDOM_next = CDOM .* depth_old ./ depth_next; %m-1
    CDOM = CDOM_next;
    
    salinity_next = (1 - ((1 - salinity ./ sal_mar) .* depth_old ./ depth_next)) .* sal_mar;
    salinity = salinity_next;
    
    d = d + 1;
end

%% then in the shelf
S = 1;  %surface area, 1 m2
MLD = 23.8;  % water depth
V = S * MLD;
SolarSpc_wl = [300:700];
CDOM_wl = [250:700];
DOC = DOC_change.DOC(22);  %mol/m3
CDOM = DOC_change{22,7:457};
DOC_amount = DOC .* V;

%%
for d = 23:days  % modelling for 2 years
    date = (d-1) + datetime('01-Jan-2019');
    month = str2num(datestr(date,'mm'));

    % pick the appropriate solar spectra
    E_month = E(:,:,month)';
    SZA_month = SZA(:,:,month)';
    
    % calculate the Kd
    CDOMforKd = CDOM(1,51:end); %only use the 300-700nm of CDOM for Kd calculation
    A = CDOMforKd + meanPABS.pabs' + W.abs';
    BB = meanbb.bb' + W.bb';
    Kd = (1 + 0.005 .* SZA_month) .* A + 4.18 .* (1 - 0.52 .* exp(-10.8 * A)) .* BB; %Lee et al.2005
    
    F_CDOM = ((1 + 0.005 .* SZA_month) .* CDOMforKd) ./ Kd;
    % generate the matrix of absorbed photon
    Q = E_month .* S .* (1 - exp(-Kd .* MLD)) .* F_CDOM; %unit: mol nm-1 s-1
    Qperday = sum(Q .* 3600);  %total amount of photon absorbed per day, unit: mol nm-1 day-1
    
    % multiply the AQY by the absorbed photon and calculate the loss of
    % DOC
    DOC_amount_loss = sum(Qperday .* AQY_DOC_spc); %mol
    DOC_amount = DOC_amount - DOC_amount_loss;
    DOC = DOC - DOC_amount_loss ./ V; % mol/m3
    DOC_change.DOC(d) = DOC;
    DOC_change.DOC_umolL(d) = DOC .* 10^6 ./ 10^3;
    DOC_change.DOC_amount(d) = DOC_amount;
    
    % the loss of CDOM
    CDOMloss = sum(Qperday .* AQY_CDOM_spc,2); % m-1 L
    CDOMloss_pervol = CDOMloss ./ (V .* 10^3); %unit: m-1
    CDOM = CDOM - CDOMloss_pervol';
    DOC_change{d,7:457} = CDOM;

end

%%
save('F:\NTU\Research\Photodegradation experiment\photodeg modelling\model results\PhotoDegModel_SG_cloud_with_mixing.mat')

