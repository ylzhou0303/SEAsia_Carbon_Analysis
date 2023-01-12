joey = readtable('C:\NTU\Research\Photodegradation experiment\Solar simulator irradiance spectrum\2020-08-18 Yongli Suntest.xlsx','sheet',1,'range','A42:AM367');

%%
figure('color','w');
for i = 2:size(joey,2)
    plot(joey{:,1},joey{:,i},'k-');hold on
end


%% only pick the spectra where the 466nm is higher than 1700 mW/m2 to calculate the mean spectrum
joey_ok = table();
joey_ok.wavelength = joey{:,1};
i = 2;
j = 2;
while i < (size(joey,2) + 1)
    if joey{82,i}>1700
        joey_ok{:,j} = joey{:,i};
        j = j+1;
    end
    i = i+1;
end

spectra_old = mean(joey_ok{:,2:end},2,'omitnan');

%%
figure;
plot(joey_ok.wavelength, spectra_old, 'k-');
xlim([300 800])



%% new lamp
yung = readtable('F:\NTU\Research\Photodegradation experiment\irradiance spectrum\2020-08-18 Yongli Suntest.xlsx','sheet',2,'range','A42:BL367');

figure('color','w');
for i = 2:size(yung,2)
    plot(yung{:,1},yung{:,i},'b-');hold on
end


%% only pick the spectra where the 466nm is higher than 1000mW/m2 to calculate the mean spectrum
yung_ok = table();
yung_ok.wavelength = yung{:,1};
i = 2;
j = 2;
while i < (size(yung,2) + 1)
    if yung{82,i}>1700
        yung_ok{:,j} = yung{:,i};
        j = j+1;
    end
    i = i+1;
end


%%
figure;
plot(yung_ok.wavelength, spectra_new, 'b-');hold on
plot(joey_ok.wavelength, spectra_old, 'k-')
xlim([300 800])
legend('new lamp','old lamp')
ylabel('irradiance (mW/m2)')


%% import the "spectra that are ok to use" of the new lamp
spc = readtable('F:\NTU\Research\Photodegradation experiment\irradiance spectrum\solar simulator spectrum.xlsx','sheet','new lamp');


%% maybe drop the fraction above 800nm
spc = spc(:,1:242);
spc.peakwv = ones(size(spc,1),1) * NaN;

%% find out the peak
for i = 1:size(spc,1)
    spc.peakwv(i) = find(spc{i,:} == max(spc{i,:}));
end

%% normalize to 760nm
nspc = spc;
nspc{1:end, 1:end} = nspc{1:end, 1:end} ./ (spc.x760 .* ones(size(nspc)));

%% plot the normalized spectra together
figure;
for i = 1:size(nspc,1)
    plot([318:2:800], nspc{i,1:end-1}, 'k-'); hold on
end



%% calculate the mean spectrum
MeanSpc = table();
MeanSpc.wl = [318:2:800]';
MeanSpc.spc = mean(spc{:,1:end-1}',2);

%% get a photon dose spectrum
plank = 6.62607004*10^-34; %Plank constant
c = 299792458;  % speed of light   m/s
avog = 6.022*10^23; % Avogadro constant
MeanSpc.NumPho = (MeanSpc.spc .* 10^-3) .* (MeanSpc.wl .* 10^-9) ./(plank * c);  %convert energy to number of photons   n m-2 s-1
MeanSpc.Einst = MeanSpc.NumPho ./ avog;       % convert from number of photons to mol of photons   mol m-2 s-1

figure;
plot(MeanSpc.wl, MeanSpc.Einst, 'k-')
ylabel('photon dose (mol m-2 s-1)')
title('spectrum of photon dose')


%% Now interpolate the value to the energy spectrum to fill the missing gap between every 2nm
MeanSpcInt = table();
MeanSpcInt.wl = [318:1:800]';
MeanSpcInt.spc = ones(size(MeanSpcInt,1),1) .* nan;

for i = 1:2:size(MeanSpcInt,1)
    MeanSpcInt.spc(i) = MeanSpc.spc (MeanSpc.wl == MeanSpcInt.wl(i));
end

for i = 2:2:size(MeanSpcInt,1)
    MeanSpcInt.spc(i) = mean(MeanSpcInt.spc([i-1 i+1]));
end

MeanSpcInt.NumPho = (MeanSpcInt.spc .* 10^-3) .* (MeanSpcInt.wl .* 10^-9) ./(plank * c);  %convert energy to number of photons   n m-2 s-1
MeanSpcInt.Einst = MeanSpcInt.NumPho ./ avog;       % convert from number of photons to mol of photons   mol m-2 s-1
 
% %% what is the total intensity between 300nm and 400nm?
% intensity300_400 = sum (MeanSpcInt.spc(1:83));
% 
% % %% because the intensity between 300nm and 400nm should be 40W/m2, we need to adjust the spectra
% % fun = @DiffCalc;
% % fun2 = @(factor)fun(MeanSpc, factor);
% % factor0 = 1;
% % [lowest,fval] = fminsearch(fun2, factor0);
% % 
% % % so the answer is that the factor should be 0.7777
% % MeanSpc.spc_cor = MeanSpc.spc .* 0.7777;
% 

%% add the spectrum below 318nm
% from the TUV model, no irradiance below 290nm, so set the irradiance at
% 290nm to 0, and interpolate the value between 290nm and 318nm
ExtSpc = table();
ExtSpc.wl = [290:317]';
ExtSpc.spc = ones(size(ExtSpc,1),1) .* nan;

k = (MeanSpcInt.spc(1) - 0)./(318-290);
b = (290*MeanSpcInt.spc(1) - 318*0) ./ (290 - 318);
ExtSpc.spc = k .* ExtSpc.wl + b;

ExtSpc.NumPho = (ExtSpc.spc .* 10^-3) .* (ExtSpc.wl .* 10^-9) ./(plank * c);  %convert energy to number of photons   n m-2 s-1
ExtSpc.Einst = ExtSpc.NumPho ./ avog;       % convert from number of photons to mol of photons   mol m-2 s-1

%% combine the two spectrum datasets to create a complete datasets
SuntestSpc = vertcat(ExtSpc, MeanSpcInt);

figure;plot(SuntestSpc.wl, SuntestSpc.Einst,'k-');
ylabel('photon dose (mol m-2 s-1)')
title('spectrum of photon dose')

%%
load('C:\NTU\Research\Photodegradation experiment\irradiance spectrum\spectrum.mat')    