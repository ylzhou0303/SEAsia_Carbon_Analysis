%% Here to model the change in DOC and CDOM for the MALUDAM experiment, to check the validity of the model

% photodegradation modelling
load ('F:\NTU\Research\Photodegradation experiment\photodeg modelling\PARs.mat');

% SolarSpc: scalar downward irrandiance just below the water surface, unit:
% mol m-2 nm-1 s-1
%%
load ('F:\NTU\Research\Photodegradation experiment\irradiance spectrum\spectrum.mat','SuntestSpc') 
load('F:\NTU\Research\Photodegradation experiment\Maludam photodegradation May-2018\AQYs.mat','MalPhbl')


%%
V = 30 .* 10^-6;  % volume of sample water, m3
S = (47/2*0.001)^2 * 3.1415926; %area of the cell, unit: m2
L = V/S; %pathlength of the sample water, m

MLD = L;  % assume that the depth is 30m
wl = [300:700];
AQY_DOC_spc = AQY_DOC(1) .* exp(-AQY_DOC(2) .* wl);
AQY_CDOM_spc = AQY_CDOM.c1 .* exp(-AQY_CDOM.c2 .* wl);
DOC = 3250 .* (10^-6);  %unit: mol/L
CDOM = MalPhbl.a{1,13:413};

days = 18;
j = 1;
DOC_change = table();
DOC_change.time = [1:days]';
DOC_change.DOC = ones(days,1)*nan;

%%
for d = 1:days  % modelling for 2 years
    date = (d-1) + datetime('01-Jan-2019');
    month = str2num(datestr(date,'mm'));
    
    % pick the appropriate solar spectra
    lim1 = 1 + ( month - 1 ) * 24;
    lim2 = lim1 + 23;
    E_surface = SolarSpc(lim1:lim2,:);
    
    
    for hr = 1:24
        E_surface_hr = SuntestSpc.Einst(11:411)';   % pick the solar spectra for this hour
        
        % calculate the Kd
%         A = CDOM.a' + meanPABS.pabs' + W.abs';
%         BB = meanbb.bb' + W.bb';
%         Kd = (1+0.005 .* E_surface_hr.SZA) .* A + 4.18 .* (1 - 0.52 .* exp(-10.8 * A)) .* BB;
%         
        Kd = CDOM;
        % calculate the fraction of light absorbed by CDOM
        F_CDOM = CDOM ./ Kd;
        
        % calculate the photon absorbed by the CDOM across the water
        % column, assuming water depth is 30m, over this hour (3600s)
        Q = E_surface_hr .* S .* F_CDOM .* ( 1 - exp( -Kd .* MLD )) .* 3600; %unit:mol nm-1
        %Q = MalPhbl.Qtotal_wl(1,11:411);
        % multiply the AQY by the absorbed photon and calculate the loss of
        % DOC and CDOM
        DOC_loss = sum(Q .* AQY_DOC_spc) ./ (V .* 1000);  %unit: mol/L
        DOC = DOC - DOC_loss;
        
        CDOM_loss_total_column = sum(Q .* AQY_CDOM_spc,2);  % the unit of the CDOM loss calculated here is  L m-1, so need to normalize to 1L
        CDOM_loss = CDOM_loss_total_column' ./ (V .* 1000); % unit: m-1
        CDOM = CDOM - CDOM_loss;
        
    end
    
    DOC_change.DOC(j) = DOC * 1000000;  %unit: umol/L
    DOC_change{j,3:403} = CDOM; %unit: m-1
    j = j + 1;
end
        


%%

% calculate the depth-resolved underwater light field
        E = table();
        E.depth = [0:30]';   % assume the depth is 30m
        E{:,2:402} = E_surface_hr{1,3:end} .* exp(-Kd .* E.depth); % calculate the 

figure;
for i = 1:31
    plot([300:700],E{i,2:end});hold on
    pause(1)
end



%% 
figure;
for i = 1:100
    plot(DOC_change{i,3:403});hold on
    pause(0.2);
end


%% compare the natural solar spectrum and the simulated spectra
load('F:\NTU\Research\Photodegradation experiment\irradiance spectrum\spectrum.mat')
figure;
hold on
plot(SuntestSpc.wl, SuntestSpc.Einst, 'r-');
for i = 1:24
    plot(wl, E_surface{i,3:end},'k-');
end

%%
figure;
plot([290:800], MalPhbl.Qtotal_wl,'k-')



%%
figure;
for i = 1:18
    hold on
    plot(DOC_change{i,3:403});
    pause(1)
end

for i = 1:7
    plot(MalPhbl.a{i,13:413},'k-','linewidth',6);
    pause(1)
end

%%
save('F:\NTU\Research\Photodegradation experiment\photodeg modelling\Maludam exp modelling.mat')