% CALCULATION BY DAY
% GENERATE A MATRIX OF MEAN ABSORBED PHOTON FOR EACH DAY, ONLY CALCULATE
% THE CHANGES IN DOC AND CDOM BY DAY, NOT BY HOUR

% MONTE CARLO SIMULATION
% photodegradation modelling
clear
location = 'Talang';  %SG or Talang
condition = 'cloud';
months = 3; % for Talang, enter 3; for SG, enter 24
season = 'winter';  % for Talang, enter spring, summer, fall or winter; for SG, just leave it empty

matname = ['PARs_',location,'_',condition,'.mat'];
load (['C:\NTU\Research\Photodegradation experiment\photodeg modelling\',matname]);

% SolarSpc: scalar downward irrandiance just below the water surface, unit:
% mol m-2 nm-1 s-1

%% set up the input data
% For this Monte Carlo simulation for error estimation, run the loop at a
% 1-month resolution

monthdays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; % number of days in each month

if isempty(season)
    month_a = 1; month_b = 24;
elseif strcmp(season,'spring')
    month_a = 2; month_b = 4;
elseif strcmp(season,'summer')
    month_a = 5; month_b = 7;
elseif strcmp(season, 'fall')
    month_a = 8; month_b = 10;
elseif strcmp(season, 'winter')
    month_a = 11; month_b = 13;
end


S = 1;  %surface area, 1 m2
MLD = MLD;  % water depth
V = S * MLD;
SolarSpc_wl = 300:700;
CDOM_wl = 250:700;

E_all = reshape(SolarSpc{:,3:end}', 401, 24, 12);  
%sunlight spectrum at each hour for each month
%one page corresponds to one month, in each page, each COLUMN corresponds to each hour, each ROW corresponds to each wavelength
SZA_all = reshape(SolarSpc.SZA, 1, 24, 12);



%% SET UP MONTE CARLO PARAMETERS
MC = struct;
MC.n = 1500;  % run times

% set up the water depth
MC.MLD = MLD + MLD_std .* randn(MC.n,1);
% remove the negative value
MC.MLD(~(MC.MLD>0)) = nan;

% set up the PABs and bb
MC.RandforP = randn(MC.n,1);  
%if you have higher PABS, you would have higher bb, so the two parameters should vary consistently
MC.PABS = meanPABS.pabs' + meanPABS.stdev' .* MC.RandforP;
figure;plot(300:700,MC.PABS);
% There are some spectra with negative values, remove them
MC.PABS(sum(MC.PABS > 0 ,2) < size(MC.PABS,2),:)=nan; 

MC.bb = meanbb.bb' + meanbb.stdev' .* MC.RandforP;
figure;plot(300:700, MC.bb);
%remove the spectra with negative values
MC.bb (sum(MC.bb > 0, 2) < size(MC.bb,2), :)=nan;


MC.A_pabs_w = MC.PABS + W.abs';
MC.BB = MC.bb + W.bb';


%% set up the AQY for DOC
AQY_DOC_spc = AQYs_DOC.Mean';
AQY_DOC_std = AQYs_DOC.Std';

MC.RandforAQY = randn(MC.n,1);
MC.AQY_DOC = AQY_DOC_spc + AQY_DOC_std .* MC.RandforAQY;
figure;plot(300:700, MC.AQY_DOC);
%remove the spectra with negative values
MC.AQY_DOC (sum(MC.AQY_DOC > 0, 2) < size(MC.AQY_DOC,2), :)=nan;
%how many spectra are removed?
sum(isnan(MC.AQY_DOC(:,1)))

%% set up the AQY for CDOM
AQY_CDOM_spc = AQY_CDOM_250_700{:,2:end};
AQY_CDOM_std = MeanCDOM_AQY_std(21:end,:); %take 250 to 700nm

MC.AQY_CDOM = ones(451,401,MC.n).*nan;
for i = 1:MC.n
    MC.AQY_CDOM(:,:,i) = AQY_CDOM_spc + AQY_CDOM_std .* MC.RandforAQY(i);
    if sum(MC.AQY_CDOM(:,:,i)<0) > 1
        MC.AQY_CDOM(:,:,i) = nan;
    end
end

figure;
for i = 1:500
    plot(300:700,MC.AQY_CDOM(451,:,i));hold on
end

% how many spectrum is removed?
sum(isnan(MC.AQY_CDOM(1,1,:)))

%% set up the data table
MC.DOC_t = ones(months,MC.n).*nan;%record changes in DOC overtime for all iterations
MC.a350_t = ones(months,MC.n).*nan; %record the changes in a350 overtime for all iterations
%MC.CDOM_t = ones(months, 451, MC.n).*nan;%record the change in CDOM spectrum over time

%% Run Monte Carlo Simulation
for n = 1:MC.n
    DOC_umolL_t = DOC_t0;  %unit: umol/L
    DOC_mol_t = DOC_umolL_t .* 10^-6 .* 10^3 .* V; %convert to mol
    CDOMspc_t = CDOM_t0.shelf';  % unit: m-1
    
    A_pas_w = MC.A_pabs_w(n,:);
    BB = MC.BB(n,:);
    AQY_DOC = MC.AQY_DOC(n,:);
    AQY_CDOM = MC.AQY_CDOM(:,:,n);
    
    w = 1;
    if isnan(A_pas_w(1,1)) || isnan(BB(1,1)) || isnan(AQY_DOC(1,1)) || isnan(AQY_CDOM(1,1))
        w = w + 1;
        continue  %exit this iteration, and move to the next iteration
    end
    


    for month = month_a:month_b  
      
        % pick the appropriate solar spectra for this month
        if month > 12
            k = month - 12;
        else
            k = month;
        end     
        E = E_all(:,:,k)';
        SZA = SZA_all(:,:,k)';
        
        % calculate the Kd
        A = A_pas_w + CDOMspc_t(51:end); %the total absorption is the sum of PABS, water absorption and CDOM. The CDOM changes over time
        Kd = (1 + 0.005 .* SZA) .* A + 4.18 .* (1 - 0.52 .* exp(-10.8 * A)) .* BB; %Lee et al.2005
        
%         % the fraction of photons absorbed by CDOM
        F_CDOM = ((1 + 0.005 .* SZA) .* CDOMspc_t(51:end)) ./ Kd;
        % generate the matrix of absorbed photon
        Q = E .* S .* (1 - exp(-Kd .* MLD)) .* F_CDOM; %unit: mol nm-1 s-1
        Q_month = sum(Q .* 3600 .* monthdays(month));  %total amount of photon absorbed per month, unit: mol nm-1 mth-1
        
        % multiply the AQY by the absorbed photon and calculate the loss of
        % DOC
        DOC_mol_loss = sum(Q_month .* AQY_DOC); %mol
        DOC_mol_t = DOC_mol_t - DOC_mol_loss; %mol
        DOC_umolL_t = (DOC_mol_t .* 10^6) ./ (V .* 10.^3);  %umol/L
        MC.DOC_t(w,n) = DOC_umolL_t;  % umol/L
        
        % the loss of CDOM
        CDOMloss = sum(Q_month .* AQY_CDOM_spc, 2); %unit: m-1 L
        CDOMloss_pervol = CDOMloss ./ (V .* 10^3); %unit: m-1
        CDOMspc_t = CDOMspc_t - CDOMloss_pervol';
        %MC.CDOM_t(w,:,n) = CDOMspc_t;
        MC.a350_t(w,n) = CDOMspc_t(101); % a350, unit:m-1
        
        w = w + 1;
    end
end

%%
figure;
plot(1:3, MC.DOC_t);
mean(MC.DOC_t(end,:),'omitnan')

figure;
plot(1:3, MC.a350_t);
mean(MC.a350_t(end,:),'omitnan')

%% calculate the standard deviation of results from all iterations
DOC_t_std = std(MC.DOC_t,0,2,'omitnan');
a350_t_std = std(MC.a350_t,0,2,'omitnan');

savename = ['MC_PhotoDegModel_',location,'_',condition,'_',season,'.mat'];
save(['C:\NTU\Research\Photodegradation experiment\photodeg modelling\model results\',savename],'DOC_t_std','a350_t_std');
% %%
% figure;
% plot([250:700],DOC_change{1:50:end,7:end});hold on
% plot([250:700],CDOM_t0.shelf,'k-')

