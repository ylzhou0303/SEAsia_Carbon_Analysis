%% do exponential fit for the Talang PABS data

%% first, examine the exponential fit for the Kusu Hantu PABS data
load('F:\NTU\Research\Photodegradation experiment\Kd\Kdmodelling.mat')

%%
KHPABSfit = table();
KHPABSfit.Date = PABS.Date;
KHPABSfit.Site = PABS.Site;
KHPABSfit.k1 = ones(size(PABS,1),1)*nan;
KHPABSfit.k2 = ones(size(PABS,1),1)*nan;


wl = [300:400];  %do the fit for 300 - 400nm
xdata = wl;
fun = @(k, xdata) (k(1).*exp(k(2).*xdata));
k0 = [1, -0.01];

for i = 1:size(KHPABSfit,1)
    ydata = PABS{i,9:109};
    k = lsqcurvefit(fun, k0, xdata, ydata);
    KHPABSfit.k1(i) = k(1);
    KHPABSfit.k2(i) = k(2);
    
    R_sq = 1-var(ydata-fun(k,xdata))/var(ydata);
    KHPABSfit.R_sq(i) = R_sq;
end


%%
% The Rsquare of the fit are mostly above 0.95, so the exponential fit
% works. We now fit the exponential curve to the Talang PABS data

load('F:\NTU\Research\Photodegradation experiment\Talang, Sarawak PABS bb\Talang_PABS_bb.mat')

%% do exponential fit for Talang PABS
TPABS_350_400 = Talang_PABS{1:51,2:6}';

TalangPABSfit = table();
TalangPABSfit.k1 = ones(5,1).*nan;
TalangPABSfit.k2 = ones(5,1).*nan;


wl = [350:400];  
xdata = wl;
fun = @(k, xdata) (k(1).*exp(k(2).*xdata));
k0 = [1, -0.01];

for i = 1:5
    ydata = TPABS_350_400(i,:);
    k = lsqcurvefit(fun, k0, xdata, ydata);
    TalangPABSfit.k1(i) = k(1);
    TalangPABSfit.k2(i) = k(2);
    
    R_sq = 1-var(ydata-fun(k,xdata))/var(ydata);
    TalangPABSfit.R_sq(i) = R_sq;
end


%%
TPABS_300_349 = TalangPABSfit.k1 .* exp(TalangPABSfit.k2 .* [300:349]);

% add the 300-350nm data to the initial data
TPABS_300_700 = horzcat(TPABS_300_349, Talang_PABS{:,2:6}');

%%
figure;
hold on
for i = 1:5
    plot([300:700],TPABS_300_700(i,:));
end
hold off

%%
figure;
for i = 1:5
    plot([350:400],TPABS_350_400(i,:),'k-');hold on
    plot([350:400],TalangPABSfit.k1(i) .* exp(TalangPABSfit.k2(i) .* [350:400]));hold off
    pause(1)
end
%% calculate the mean and std
clear Talang_MeanPabs
Talang_MeanPabs = table();
Talang_MeanPabs.wl = [300:700]';

joey = mean(TPABS_300_700,1);
Talang_MeanPabs.PABS = joey';

yung = std(TPABS_300_700,0,1);
Talang_MeanPabs.stdev = yung';

%%
figure;
hold on
plot(Talang_MeanPabs.wl, Talang_MeanPabs.PABS,'k-');
plot(Talang_MeanPabs.wl, Talang_MeanPabs.PABS + Talang_MeanPabs.stdev, 'b-');
plot(Talang_MeanPabs.wl, Talang_MeanPabs.PABS - Talang_MeanPabs.stdev, 'b-');


%%
save('F:\NTU\Research\Photodegradation experiment\Talang, Sarawak PABS bb\Talang_PABSfit_bb.mat')
