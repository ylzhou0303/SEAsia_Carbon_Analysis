%%% This file calculates the mean and std for AQY for CDOM and DOC

%% CDOM AQY
%% %% AQY for CDOM
load('F:\NTU\Research\Photodegradation experiment\Experimental Data\Maludam photodegradation May-2018\AQYs.mat')
CDOM_AQY{CDOM_AQY.aloss<0, 3:4} = NaN; %remove the data point with negative CDOM loss (i.e., CDOM increase after photodeg exp)
Mal_CDOM_AQY = CDOM_AQY;
clear CDOM_AQY


load('F:\NTU\Research\Photodegradation experiment\Experimental Data\Maludam Mix with seawter\Mal_mix_AQYs.mat')
CDOM_AQY{CDOM_AQY.aloss<0, 3:4} = NaN;
MalMix_CDOM_AQY = CDOM_AQY;
clear CDOM_AQY

load('F:\NTU\Research\Photodegradation experiment\Experimental Data\Kusu phbl\Kusu Jul 2020\KusuPhblAQY.mat')
CDOM_AQY{CDOM_AQY.aloss<0, 3:4} = NaN;
Kusu_CDOM_AQY = CDOM_AQY;
clear CDOM_AQY

load('F:\NTU\Research\Photodegradation experiment\Experimental Data\20201105 optical filter\Mal_OF_AQYs.mat')
% use the AQY calculated from all optical filter treatments
CDOM_AQY_OF{272:end, 2:4} = NaN;  % set AQY above 500nm NaN, because those wavelengths have lots of noise
MalOF_CDOM_AQY = CDOM_AQY_OF;
clear CDOM_AQY


clearvars -except Mal_CDOM_AQY MalOF_CDOM_AQY MalMix_CDOM_AQY Kusu_CDOM_AQY

%%
figure;
mesh([230:700],[300:700],MeanCDOM_AQY')

%%
figure;
plot(squeeze(CDOM_AQY_3D(121,:,1:4)))
legend('Mal','MalMix','Kusu','Mal_OF')
set(gca, 'yscale','log')

%%
CDOM_AQY_3D = ones(471,401, 4) .* nan;
CDOM_AQY_3D(:,:,1) = Mal_CDOM_AQY.cd_1 .* exp( - Mal_CDOM_AQY.cd_2 .* [300:700]);
CDOM_AQY_3D(:,:,2) = MalMix_CDOM_AQY.cd_1 .* exp( - MalMix_CDOM_AQY.cd_2 .* [300:700]);
CDOM_AQY_3D(:,:,3) = Kusu_CDOM_AQY.cd_1 .* exp( - Kusu_CDOM_AQY.cd_2 .* [300:700]);
CDOM_AQY_3D(:,:,4) = MalOF_CDOM_AQY.cd_1 .* exp( - MalOF_CDOM_AQY.cd_2 .* [300:700]);


%% the mean spectrum of CDOM AQY, dimension: 471*401, i.e. 230:700 * 300:700
MeanCDOM_AQY = mean(CDOM_AQY_3D,3,'omitnan');
MeanCDOM_AQY_std = std(CDOM_AQY_3D,0,3,'omitnan');

%%
save('F:\NTU\Research\Photodegradation experiment\photodeg modelling\MeanCDOMAQY.mat','CDOM_AQY_3D','MeanCDOM_AQY','MeanCDOM_AQY_std')



%% DOC AQY
AQYs_DOC = table();

load('F:\NTU\Research\Photodegradation experiment\Experimental Data\Maludam photodegradation May-2018\AQYs.mat')
wl = [300:700]';
AQYs_DOC.Exp1 = cd_opt(1) .* exp(-cd_opt(2) .* wl);
clearvars -except AQYs_DOC
 
load('F:\NTU\Research\Photodegradation experiment\Experimental Data\Maludam Mix with seawter\Mal_mix_AQYs.mat')
wl = [300:700]';
AQYs_DOC.Exp2 = cd_opt(1) .* exp(-cd_opt(2) .*wl);
clearvars -except AQYs_DOC 

load('F:\NTU\Research\Photodegradation experiment\Experimental Data\Kusu phbl\Kusu Jul 2020\KusuPhblAQY.mat')
wl = [300:700]';
AQYs_DOC.Exp3 = cd_opt(1) .* exp(-cd_opt(2) .* wl);
clearvars -except AQYs_DOC 

load('F:\NTU\Research\Photodegradation experiment\Experimental Data\20201105 optical filter\Mal_OF_AQYs.mat')
wl = [300:700]';
AQYs_DOC.Exp4 = c2(1) .* exp(-c2(2) .* wl);
clearvars -except AQYs_DOC 
% for Exp4, here we use the optimized AQY calculated using the multiple
% optical filter treatments following Powers et al. 2017

%% calculate the mean and std
AQYs_DOC.Mean = mean(AQYs_DOC{:,1:4},2);
AQYs_DOC.Std = std(AQYs_DOC{:,1:4},0,2);

load('F:\NTU\Research\Photodegradation experiment\photodeg modelling\AQYs_DOC.mat','AQYs_DOC');


%% plot the AQYs
figure;
wl = [300:700];

plot(wl, AQYs_DOC.Exp1,'r--');hold on
plot(wl, AQYs_DOC.Exp2,'g--');
plot(wl, AQYs_DOC.Exp3,'b--');
plot(wl, AQYs_DOC.Exp4,'y--')
plot(wl, AQYs_DOC.Mean,'k-');
plot(wl, AQYs_DOC.Mean + AQYs_DOC.Std,'b-');
plot(wl, AQYs_DOC.Mean - AQYs_DOC.Std,'b-');

ylim([0 0.001])
set(gca,'ytick',[0 2 4 6 8 10].*10^-4)
%%
for i = 2:5
    c = AQYs_DOC{1,i};
    d = AQYs_DOC{2,i};
    plot(wl, c.*exp(-d.*wl),'--');
end

