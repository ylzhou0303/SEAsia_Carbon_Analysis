%% In this file, I use a power law to fit the back scattering data and extrapolate to below 400nm, then reconstruct the Kd and compare it with the measured Kd

% import the back scattering data
bb = readtable('F:\NTU\Research\Photodegradation experiment\Kd\Final dataset Matlab.xlsx','sheet','Particulate backscattering','range','A3:N39');

%% Fit the power law to each bb spectrum
bb.k1 = ones(size(bb,1),1)*nan;
bb.k2 = ones(size(bb,1),1)*nan;

wl = [412, 440, 488, 510, 532, 595, 650, 676, 715];
xdata = wl;
fun = @(k,xdata) k(1).*(xdata.^k(2));  %define a power law function for the curvefitting
k0 = [1,-0.2];  %

for i = 1:size(bb,1)
    ydata = bb{i,6:14};
    [k] = lsqcurvefit(fun, k0, xdata, ydata);
    bb.k1(i) = k(1);
    bb.k2(i) = k(2);
    
    R_sq = 1-var(ydata-fun(k,xdata))/var(ydata);
    bb.R_sq(i) = R_sq;
end


%% plot the fitted curve
figure;
for i = 1:size(bb,1)
    plot(wl,bb{i,6:14},'ko');hold on
    plot(wl,bb.k1(i).*wl.^bb.k2(i),'b-');hold off
    pause(2)
end

%% look good, now do the interpolation and extrapolation
wl = [300:1:700];
bbfit = table();
bbfit.Date = bb.Date;
bbfit.Site = bb.Site;
bbfit.Station = bb.Station;

bbfit{:,4:404} = bb.k1 .* wl .^ bb.k2;  %these are the interpolated and extrapolated bb values. 300 - 900nm

% the R_sq of the 2020 data is poor, so do not use the 2020 bb data
bbfit = bbfit(1:24,:);
%%
save('F:\NTU\Research\Photodegradation experiment\Kd\bbcurvefitting.mat')