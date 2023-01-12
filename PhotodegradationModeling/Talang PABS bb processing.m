%% compile the Talang PABS and bb

Talang_PABS = readtable('F:\NTU\Research\Photodegradation experiment\Talang, Sarawak PABS bb\Talang_PABS_bb.xlsx','sheet','raw data');

%% identify problematic spectra
figure;
wl = round(Talang_PABS.wl);

hold on
for i = 3:2:size(Talang_PABS,2)
    plot(wl,Talang_PABS{:,i});
    pause(2)
end


%% identify problematic spectra
figure;
wl = round(Talang_PABS.wl);

hold on
for i = 4:2:size(Talang_PABS,2)
    plot(wl,Talang_PABS{:,i});
    pause(2)
end

%no weired spectra identified

%% BASED ON THE LOCATION OF THE STATIONS, WE TAKE ST4,10,11,12,16 for our calculation. These stations are at the coast.
Talang_PABS_raw = readtable('F:\NTU\Research\Photodegradation experiment\Talang, Sarawak PABS bb\Talang_PABS_bb.xlsx','sheet','raw data');


%%
Talang_ap = Talang_PABS_raw(:,[1 2:2:10]);  %phytoplankton absorption
Talang_ad = Talang_PABS_raw(:,[1 3:2:11]);  %detritus absorption

Talang_ap.wl = round(Talang_ap.wl);
Talang_ad.wl = round(Talang_ad.wl);

%% add the phytonplankton absorption and detritus absorption together as the PABS
Talang_PABS = table();
Talang_PABS.wl = [350:700]';

% because the wavelength of the raw data does not match our purpose, we
% need to process it a bit.
for i = 1:size(Talang_PABS,1)
    k = find(Talang_ap.wl == Talang_PABS.wl(i));
    
    ap = mean(Talang_ap{k,2:end},1);
    ad = mean(Talang_ad{k,2:end},1);
    pabs = ap + ad;
    
    Talang_PABS{i,2:6} = pabs;
end


%% calculate the mean Talang PABS and stdev for each wavelength
Talang_MeanPabs = table();
Talang_MeanPabs.wl = [350:700]';
Talang_MeanPabs.PABS = mean(Talang_PABS{:,2:6},2);
Talang_MeanPabs.stdev = std(Talang_PABS{:,2:6},0,2);


%%
figure;
hold on
for i = 2:size(Talang_PABS,2)
    plot(Talang_PABS{:,i});
    pause(1);
end
hold off
%%
figure;hold on
plot(Talang_MeanPabs.wl, Talang_MeanPabs.PABS,'k-');
plot(Talang_MeanPabs.wl, Talang_MeanPabs.PABS + Talang_MeanPabs.stdev, 'b-');
plot(Talang_MeanPabs.wl, Talang_MeanPabs.PABS - Talang_MeanPabs.stdev, 'b-');


%% process the bb data
Talang_bb = readtable('F:\NTU\Research\Photodegradation experiment\Talang, Sarawak PABS bb\Talang_PABS_bb.xlsx','sheet','bb');

Talang_bb = Talang_bb(1:6:65,:);


%%
figure;hold on
for i = 1:11
    plot(wl, Talang_bb{i,2:10});
    pause(1)
end

%% step 1: fit the power law curve to the data
Talang_bb.k1 = ones(size(Talang_bb,1),1)*nan;
Talang_bb.k2 = ones(size(Talang_bb,1),1)*nan;

wl = [412, 440, 488, 510, 532, 595, 650, 676, 715];
xdata = wl;
fun = @(k,xdata) k(1).*(xdata.^k(2));  %define a power law function for the curvefitting
k0 = [1,-0.2];  %

for i = 1:size(Talang_bb,1)
    ydata = Talang_bb{i,2:10};
    k = lsqcurvefit(fun, k0, xdata, ydata);
    Talang_bb.k1(i) = k(1);
    Talang_bb.k2(i) = k(2);
    
    R_sq = 1-var(ydata-fun(k,xdata))/var(ydata);
    Talang_bb.R_sq(i) = R_sq;
    clear k
end


%%
figure;
wl = [300:700];
plot(wl,Talang_bb.k1(1).*(wl .^ Talang_bb.k2(1)),'b-');hold on
plot(wl,6.8399.*(wl .^ (-0.8726)),'r-');
plot([412 440 488 510 532 595 650 676 715], Talang_bb{1,2:10},'ko')



%% create the full spectra and calculate the mean and stdev
Talang_bbfit = table();
Talang_bbfit.Stn_no = Talang_bb.Stn_no;
Talang_bbfit{:,2:402} = Talang_bb.k1 .* ([300:700] .^ Talang_bb.k2);


%% check the model fit
figure;
for i = 1:11
    plot([300:700], Talang_bbfit{i,2:402},'k-');hold on
    plot([412, 440, 488, 510, 532, 595, 650, 676, 715], Talang_bb{i,2:10},'ko');
    hold off
    pause(2)
end


%% looks good, calculate the mean bb and stdev
Talang_Meanbb = table();
Talang_Meanbb.wl = [300:700]';
Talang_Meanbb.bb = mean(Talang_bbfit{:,2:end}',2);
Talang_Meanbb.stdev = std(Talang_bbfit{:,2:end}',0,2);

%%
figure;
plot(Talang_Meanbb.wl, Talang_Meanbb.bb, 'k-');hold on
plot(Talang_Meanbb.wl, Talang_Meanbb.bb + Talang_Meanbb.stdev, 'b-');
plot(Talang_Meanbb.wl, Talang_Meanbb.bb - Talang_Meanbb.stdev, 'b-');

%%

save('F:\NTU\Research\Photodegradation experiment\Talang, Sarawak PABS bb\Talang_PABS_bb.mat')