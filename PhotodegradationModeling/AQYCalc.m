% here I need to work out the AQY for the Maludam photodegradation
% experiment

% load the photodegradation data
MalPhblData = readtable ('F:\NTU\Research\Photodegradation experiment\Experimental Data\Maludam photodegradation May-2018\Maludam phbl experiment Jun 2018 DOC CDOM data.xlsx','sheet','light');

%load the spectrum of the solar simulator the unit of the photo flux is mol m-2 s-1 nm-1
load ('F:\NTU\Research\Photodegradation experiment\Solar simulator irradiance spectrum\solar simulator spectrum.mat','SuntestSpc') 
load('F:\NTU\Research\Photodegradation experiment\Kd\Kdmodelling.mat','W'); %water absorption

%% Tidy up the absorption coefficients
MalPhbl = struct('import_data',MalPhblData);

%calculate the mean CDOM spectra for the starting and end point
MalPhbl.a_raw = MalPhbl.import_data;  
MalPhbl.a_raw(13:15,:) = [];
MalPhbl.a_raw{1,:} = mean(MalPhbl.import_data{1:3,:},1,'omitnan'); %calculate the mean for T0

MalPhbl.a_raw{2:11,:} = MalPhbl.import_data{4:13,:};
MalPhbl.a_raw{12,:} = mean(MalPhbl.import_data{14:15,:},1,'omitnan'); %calculate the mean for T_end

%% extract the timepoint, DOC, and spectra between 290nm and 700nm and 
%define which part of the data you want to use, here I use 0-431 hours, i.e. row 1-7
%because during the rest of the experiment, flocs were observed
MalPhbl.a = MalPhbl.a_raw(1:7,[2:3,73:483]);  
MalPhbl.a{:,3:end} = MalPhbl.a{:,3:end} .* (MalPhbl.a{:,3:end} > 0);  %set all the negative absorption coefficients to 0


%% calculate the mean absorption for each time section
MalPhbl.mean_a = MalPhbl.a;
MalPhbl.duration = MalPhbl.mean_a.time(2:end) - MalPhbl.mean_a.time(1:end-1);

for i = 2:size(MalPhbl.mean_a,1)
    MalPhbl.mean_a{i,2:end} = mean(MalPhbl.a{i-1:i,2:end});
end
MalPhbl.mean_a(1,:) = [];


%% calculate the AQY using the data from 0 - 431 hour, because after that, flocs were observed until the experiment was terminated
% calculate the photons absorbed during each time section, e.g. 0-72hour, individually

% calculate the amount of photons that's absorbed by the sample (throughout the entire sample cell) during each time section
V = 30 .* 10^-6;  % volume of sample water, m3
S = (47/2*0.001)^2 * 3.1415926; %area of the cell, unit: m2
PL = V/S; %pathlength of the sample water, m

CDOM_Napierian_a = MalPhbl.mean_a{:,3:end};
Kd = W.abs(91:501)' + W.bb(91:501)' + CDOM_Napierian_a;  % the Kd is the sum of water absorption, backscattering and CDOM
F_CDOM = (1 - exp(-Kd .* PL)) .*  (CDOM_Napierian_a ./ Kd);   % the fraction of irradiance absorbed by CDOM across the entire water column

%%
MalPhbl.Q = MalPhbl.mean_a; 
MalPhbl.Q{:,3:end} = nan;
MalPhbl.Q{:,3:end} = SuntestSpc.Einst(1:411)' .* 0.95 .* S .* F_CDOM;  

%unit: mol photons s-1
% the transmission of light through the quartz material is accounted for by
% a multiplying a factor of 0.95. The material of our quartz cells is
% Spectrosil® Quartz. The transmission of light is 95%. (https://www.starna.com/cells/cell-specifications)

% in this step, the light spectrum is multiplied by the surface area of the
% cuvette, to get the amount of photons that fall onto the surface of the
% cuvette. It is then multiplied the fraction of the absorbed incident light.
% The fraction is calculated based on this equation:
% A = -lg (Pt - Pi), where A is the absorbance (m-1) and Pt is the
% transmitted light and Pi is the incident light, the fraction of absorbed
% light would then be: 1 - (Pt/Pi) = 1 - 10^(-A)

% multiply the time and calculate the total amount of photons absorbed
MalPhbl.Qtotal_wl = sum(MalPhbl.Q{:,3:end} .* MalPhbl.duration .* 3600, 'omitnan'); % spectrally resolved abosrbed photons
MalPhbl.Qtotal = sum(MalPhbl.Qtotal_wl,'omitnan');  % absorbed photons integrated over the spectral range

%% now calculate the AQY
% spectrum-integrated AQY
DOCloss = MalPhbl.a.DOC(1) - MalPhbl.a.DOC(end);   % umol/L
DOCloss = (DOCloss .* 10^-6 .* 10^3) .* V;      % mol
fai = DOCloss ./ MalPhbl.Qtotal;   %mol C (mol photon -1)


%% calculate the AQY using the data from different time sections
BroadBandAQY = table();
BroadBandAQY.time_section = {'0-72','72-144','144-216','216-288','288-358','358-431'}';

for k = 1:size(MalPhbl.Q,1)
   BroadBandAQY.DOCloss(k) = (MalPhbl.a.DOC(k) - MalPhbl.a.DOC(k+1)) * 10^-6 .* 10^3 .* V;
   BroadBandAQY.Q(k) = sum(MalPhbl.Q{k,3:end} .* MalPhbl.duration(k) .* 3600);
end
BroadBandAQY.AQY = BroadBandAQY.DOCloss ./ BroadBandAQY.Q .* 10^6;




%% now calculate the spectrally resolved AQY using fminsearch function
c0 = [1, 0.04];
MalPhbl.Qtotal_spc = MalPhbl.Qtotal_wl;
fun = @(c)AQYcalc(c, MalPhbl.Qtotal_spc, DOCloss);
[c,fval] = fminsearch(fun, c0);


% plot the spcetral AQY
AQYspc = c(1) .* exp(-c(2).* [290:700]);
figure;plot([290:700],AQYspc,'k-');hold on
title('spectral AQY')
xlabel('wavelength(nm)')
ylabel('AQY (mol C (mol photon)-1 nm(-1))')


%% now calculate the spectrally resolved AQY using fminsearch function (use a big c starting value)
c0 = [100, 0.04];
MalPhbl.Qtotal_spc = MalPhbl.Qtotal_wl;
fun = @(c)AQYcalc(c, MalPhbl.Qtotal_spc, DOCloss);
[c_big,fval_big] = fminsearch(fun, c0);


%% now calculate the spectrally resolved AQY using fminsearch function (use a small c starting value)
c0 = [0.01, 0.04];
MalPhbl.Qtotal_spc = MalPhbl.Qtotal_wl;
fun = @(c)AQYcalc(c, MalPhbl.Qtotal_spc, DOCloss);
[c_small,fval_small] = fminsearch(fun, c0);

%%
% plot the spcetral AQY
AQYspc = c(1) .* exp(-c(2).* [290:700]);
AQYspc_small = c_small(1) .* exp(-c_small(2).* [290:700]);
AQYspc_big = c_big(1) .* exp(-c_big(2).* [290:700]);

figure;
plot([290:700],AQYspc,'k-');hold on
plot([290:700],AQYspc_big,'r-');
plot([290:700],AQYspc_small,'b-');

legend('c0=1','c0=100','c0=0.01')
title('spectral AQY')
xlabel('wavelength(nm)')
ylabel('AQY (mol C (mol photon)-1 nm(-1))')

%% because the optimized c and d values vary with the initial c,d,
%we need a monte carlo simulation to find out the combination of c and d
%with the lowest residual

cd = zeros(1,2);
cd_opt = zeros(1,2);
wl = 290:1:700;
sl = 300:1:500;
c0 = 1.0;
d0l = 0.011;
d0h = 0.084;
montestps = int16(1000); % number of iterations
loc = 95;
lcnt = int16(length(wl));
fvals = zeros(montestps, lcnt + 4);
optarg = optimset('Display','off');
f_opt = zeros(1, lcnt);
cnf_int = zeros(2, lcnt);
loc_eps = int16((((100.0 - loc)/100.0)*montestps)/2.0);

f_var = zeros(2, lcnt);
f_var(1, 1:lcnt) = realmin;
f_var(2, 1:lcnt) = realmax;

%% iteration of optimization
for n = 1:montestps
    cd_opt = [-1,-1]; 
    exitflag = 0;
   
    while (cd_opt(1)<0 || cd_opt(2)<0 || exitflag~= 1)  %when the optimized c or d is negative, the optimization is repeated
        cd(1) = c0;
        cd(2) = d0l + ((d0h - d0l) * rand);
        [cd_opt, fval, exitflag, output] = fminsearch(@(x)AQYcalc(x,MalPhbl.Qtotal_wl,DOCloss), cd, optarg); 
    end
    
    [ress, pmzest, fitl] = AQYcalc(cd_opt, MalPhbl.Qtotal_wl, DOCloss);
    fvals(n,1:lcnt) = fitl;   % the 1:lcnt columns are the modelled fi at each wavelength
    fvals(n,lcnt+1) = cd_opt(1);
    fvals(n,lcnt+2) = cd_opt(2);
    fvals(n,lcnt+4) = ress;
end


% pick the optimization with the lowest ress
k = find(fvals(:,lcnt+4) == min(fvals(:,lcnt+4)));
cd_opt = zeros(1,2);
cd_opt(1) = fvals(k,lcnt+1);
cd_opt(2) = fvals(k,lcnt+2);


%% calculate uncertainty
[fvals_sort, sort_index] = sort(fvals,'descend');
cnf_int(1,1:lcnt) = fvals_sort(loc_eps,1:lcnt);    % below the first 5%
cnf_int(2,1:lcnt) = fvals_sort((montestps - loc_eps),1:lcnt);    %above the last 5%
% this would constrain the output between 5% and 95%

%% plot the AQY with uncertainty
figure('color','w')
x = [300:700];
y = cd_opt(1) .* exp(-cd_opt(2).*x);
semilogy(x,y,'k-');%the optimized AQY
hold on
semilogy(x,cnf_int(1,11:end),'k--');
semilogy(x,cnf_int(2,11:end),'k--')
ylim([10^-9 10^-3])

%% plot all the attempts
figure;
semilogy(x,fvals(:,11:lcnt),'k-');hold on
plot(x,y,'r-')

%% compare the optimized AQY and the preliminary AQY
figure;
plot([290:700], c(1) .* exp(-c(2).* [290:700]),'k-');hold on
plot([290:700], cd_opt(1) .* exp(-cd_opt(2) .* [290:700]),'r-')

%% calculate the AQY of CDOM
%calculate the AQY for the absorption coefficient for each wavelength

CDOM_AQY = table();
CDOM_AQY.CDOMwl = [230:700]';
CDOM_AQY.aloss = (MalPhbl.a_raw{1,13:483}' - MalPhbl.a_raw{7, 13:483}') .* 30 ./ 1000;  
% the loss of absorption coefficient multiplied by the volume, unit:L m-1
% only use the data betweeen hour 0 and hour 431

%% using the fminsearch function to calculate the spectral AQY for each CDOM wavelength
CDOM_AQY.cd_1 = ones(length(CDOM_AQY.CDOMwl),1) .* nan;
CDOM_AQY.cd_2 = ones(length(CDOM_AQY.CDOMwl),1) .* nan;

for k = 1:size(CDOM_AQY,1)
    cd0 = [100, 0.004];
    fun = @(cd)AQYcalc(cd, MalPhbl.Qtotal_wl, CDOM_AQY.aloss(k));
    [cd,fval] = fminsearch(fun, cd0);
    CDOM_AQY.cd_1(k) = cd(1);
    CDOM_AQY.cd_2(k) = cd(2);
end


%%
figure;
plot(cd(1).*exp(-cd(2).*[290:700]))

%% test, reconstruct the loss of CDOM using the modelled AQY
CDOM_AQY.aloss_calc = ones(size(CDOM_AQY,1),1) .* nan;
wl = [290:700]';

for i = 1:size(CDOM_AQY,1)
    spec_aqy = CDOM_AQY.cd_1(i) .* exp( - CDOM_AQY.cd_2(i) .* wl);
    joey = spec_aqy .* MalPhbl.Qtotal_spc';
    CDOM_AQY.aloss_calc(i) = sum(joey);
end

%%
% the reconstructed loss of CDOM between 230nm and 700nm look great, which
% means that the modelled AQY can be used for modelling the photobleaching
% of this part of the CDOM spectrum. However, the data above 700nm are
% problematic, which is mainly due to the very low amount of abosrbed
% photons. In addition, because CDOM should be zero beyond 700nm, so the
% wavelength range of CDOM for photodeg modelling is 290-700nm.

%% make some plots
% spectrum of absorbed photons
figure;plot([290:800],MalPhbl.Qtotal_spc,'k-');
xlabel('wavelength(nm)')
ylabel('photon absorbed (mol m-2 nm-1)')

%% test the difference between with and without 290-300nm
AA = sum(AQYspc .* MalPhbl.Qtotal_wl);
BB = sum(AQYspc(1,10:end) .* MalPhbl.Qtotal_wl(1,10:end));

(BB-AA)/AA
% = 0.9835


%% use the optimized AQY to reproduce DOC and CDOM loss, and compare with the measured value
wl = [290:700];
AQY_spc = cd_opt(1) .* exp (-cd_opt(2) .* wl);
DOC_pred = sum(AQY_spc .* MalPhbl.Qtotal_wl);

reldiff = ((DOC_pred - DOCloss)./DOCloss);

%% calculate the predicted DOC for each time section
Pred = table();
Pred.time_section = {'0-72','72-144','144-216','216-288','288-358','358-431'}';
Pred.DOCloss_meas = BroadBandAQY.DOCloss;

Pred.DOCloss_pred = sum(MalPhbl.Q{:,3:end} .* MalPhbl.duration .* 3600 .* AQY_spc, 2);
Pred.Q =  sum(MalPhbl.Q{:,3:end} .* MalPhbl.duration .* 3600, 2);
nRMSE = (sum((Pred.DOCloss_pred - Pred.DOCloss_meas) .^ 2) ./ 6) .^ 0.5 ./ (mean(Pred.DOCloss_meas));
Pred.reldiff = ((Pred.DOCloss_pred - Pred.DOCloss_meas) ./ Pred.DOCloss_meas) ;

figure;
plot(Pred.DOCloss_meas, Pred.DOCloss_pred, 'ko');hold on
plot(DOCloss, DOC_pred, 'ro')
plot([10^-6 7*10^-5], [10^-6 7*10^-5],'k-')



%% 
save('C:\NTU\Research\Photodegradation experiment\Experimental Data\Maludam photodegradation May-2018\AQYs.mat')



%% nested functions 
function [ress,DIC_est,AQY] = AQYcalc (x, Qtotal_spc, DOCloss)
    %calculate the DIC production
    AQY = x(1) .* exp( - (x(2) .* [290:1:700]));   % spectral AQY    
    DIC_est = sum(AQY .* Qtotal_spc,'omitnan');
    diff = DOCloss - DIC_est;
    relative_diff = diff ./ DOCloss;
    ress = relative_diff.^2;
end