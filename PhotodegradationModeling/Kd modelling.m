%% This file is to reconstruct the Kd based on the PABS, CDOM, extended bb, bb of water and absorption of water
% import the extended bb
load('F:\NTU\Research\Photodegradation experiment\Kd\bbcurvefitting.mat','bbfit')

CDOM = readtable('F:\NTU\Research\Photodegradation experiment\Kd\Final dataset Matlab.xlsx','sheet','CDOM, S, T, Chl-a','range','A5:ZH134');
PABS_raw = readtable('F:\NTU\Research\Photodegradation experiment\Kd\Final dataset Matlab.xlsx','sheet','Particulate absorption','range','A4:WK120');
Kd_meas_raw = readtable('F:\NTU\Research\Photodegradation experiment\Kd\1 Kd Singapore all spectra for MATLAB.xlsx','range','A5:LQ67');
W = readtable('F:\NTU\Research\Photodegradation experiment\Kd\water a and b MATLAB.xlsx','sheet','MATLAB');

%%
Kd_meas = Kd_meas_raw(1:2:end,:); %only keep the actual spectra, and only from 318nm to 900nm
PABS = PABS_raw(strcmp(PABS_raw.Type,'total'),:);  %only keep the total PABS, remove the detrital, from 318-900nm
W.bb = W.b ./ 2;   %calculate the water backscattering from the total scattering

%% uniform the date format
CDOM.Date = cellstr(datestr(CDOM.Date,'dd-mmm-yyyy'));
PABS.Date = cellstr(datestr(PABS.Date,'dd-mmm-yyyy'));
Kd_meas.datetime = cellstr(datestr(Kd_meas.datetime, 'dd-mmm-yyyy'));
bbfit.Date = cellstr(datestr(bbfit.Date, 'dd-mmm-yyyy'));


%%
Kd_mod = table();
Kd_mod.Date = Kd_meas.datetime;
Kd_mod.Site = Kd_meas.site;
Kd_mod.Station = Kd_meas.station;

for i = 1:size(Kd_meas)
    cdom = ones(1,292)*nan;  % to match the resolution of the Kd_meas, which is 318nm to 900nm with interval of 2nm.
    pabs = ones(1,292)*nan;
    bb = ones(1,292)*nan;
    a = ones(1,292)*nan;
    b = ones(1,292)*nan;

    
    for k = 1:size(CDOM,1)
        if strcmp(char(CDOM.Date(k)), char(Kd_meas.datetime(i))) && strcmp(CDOM.Site(k),Kd_meas.site(i)) 
            cdom = CDOM{k,102:2:684};
        end
    end
    
    for k = 1:size(PABS,1)
        if strcmp(char(PABS.Date(k)), char(Kd_meas.datetime(i))) && strcmp(PABS.Site(k), Kd_meas.site(i)) 
            pabs = PABS{k,27:2:end};
        end
    end
    
    for k = 1:size(bbfit,1)
        if strcmp(char(bbfit.Date(k)), char(Kd_meas.datetime(i))) && strcmp(bbfit.Site(k), Kd_meas.site(i)) 
            bb = bbfit{k,22:2:end};
        end
    end
    
    % now reconstruct the Kd
    a = cdom + pabs + W.abs(119:2:701)';
    b =  bb + W.bb(119:2:701)';
    
    Kd_mod{i,4:295} = (1+0.005 .* Kd_meas.SZA(i)) .* a + 4.18 .* (1 - 0.52 .* exp(-10.8 * a)) .* b;  %Lee et al. 2005 (eq. 11)
    
end



%% OK. Let's compare the model results with the measured data
figure;
for i = 1:size(Kd_meas,1)
    if ~isnan(Kd_mod{i,4})
        plot([318:2:950], Kd_meas{i,13:end},'r-');hold on
        plot([318:2:900], Kd_mod{i,4:end},'k-');hold off
        legend('measured','modelled')
        xlabel('nm')
        ylabel('Kd')
        pause(2)
    end
end
    
%%
figure;plot(W.abs)
figure;
for i = 1:size(PABS,1)
    plot([300:900],PABS{i,9:end});
    pause(0.5)
end

%%
save('F:\NTU\Research\Photodegradation experiment\Kd\Kdmodelling.mat')