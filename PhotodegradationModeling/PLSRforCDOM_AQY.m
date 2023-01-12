load spectra
X = NIR;
y = octane;

[XL,yl,XS,YS,beta,PCTVAR] = plsregress(X,y,10);

y_pred = sum(beta' .* horzcat(1, NIR(1,:)));


%%
% In this file, I calculate the AQY for CDOM photobleaching following Zhu
% et al. 2020, using PLSR + Gradient descent algorithm

%% load data of cut-off filter experiment
load('C:\NTU\Research\Photodegradation experiment\Experimental Data\20201105 optical filter\Mal_OF_AQYs.mat','Data_raw','Data','W','SuntestSpc');
Data_raw = Data_raw([1:6 8:15],:);

% Irradiance and CDOM absorption and absorbed photons:290 - 700nm
% CDOMloss: 230 - 700nm

%% Calculate the CDOM loss, which is the response variable, y
CDOMloss = Data_raw{1, 15:485} - Data_raw{2:12, 15:485};  %take CDOM from 230nm to 700nm

%convert the data from m-1 to m2 (m3 * m-1)
V = 30 .* 10^-6;  %sample volume, 30 ml, unit: m3
CDOMloss = CDOMloss .* V;

%% Calculate the absorbed photons, which is the independent variable, X

% the mean CDOM absorption during experiment, used for calculating photons
% absorbed during experiment
CDOMmean = zeros(11,411) * nan;
n = 2;
for k = 1:11
    CDOMmean(k,:) = mean(vertcat(Data_raw{1, 75:485}, Data_raw{n, 75:485}) ,1);
    n = n + 1;
end

%% calculate the light spectrum for each filter treatment
LightSpc = zeros(11,411) * nan;
Cut_off = Data_raw.filter(2:12);

for k = 1:11
    %  create an array of transmittance dependent on the filter wl, the transmittance above the filter wl is 95%
    % the optical filter can filter out the light below its specified wl,
    % and allow the light above its specified wl to pass
    T = (SuntestSpc.wl(1:411) > Cut_off(k)) .* 0.95;
    
    % then multiply the solar simulator spectrum by the T
    LightSpc(k, :) = (SuntestSpc.Einst(1:411) .* T)';
end


%% calculate the number of photons absorbed
Q = zeros(11,401) * nan;
V = 30 .* 10^-6;  % volume of sample water, m3
S = (47/2*0.001)^2 * 3.1415926; %area of the cell, unit: m2
PL = V/S; %pathlength of the sample water, m
time = 144 * 3600; % duration of irradiation, unit: s

Kd = W.abs(91:501)' + W.bb(91:501)' + CDOMmean;  % the Kd is the sum of water absorption, backscattering and CDOM
F_CDOM = (1 - exp( - Kd .* PL )) .* CDOMmean ./ Kd; %the fraction of CDOM in the total Kd
Q = LightSpc .* 0.95 .* S .* time .* F_CDOM;
% spectrally resolved absorbed photon (absorbed by CDOM), unit: mol nm-1
%0.95 is the transmittance of the glass of the cuvette


%% Now run the PSLR for CDOM absorption at each wavelength
Beta = zeros(412,431) .* nan;

for i = 1:471
    X = Q;
    y = CDOMloss(:,i);

    [XL,yl,XS,YS,beta,PCTVAR] = plsregress(X,y,5);
    Beta(:,i) = beta;
end

% first row is intercept
% In the Beta table, when looking at a specific column, that is a specific 
% wavelength of CDOM absorption, from top to bottom is the coefficient of absorbed photons at each wavelength of irradiance  


%% Plot the percent of variance explained in the response variable (PCTVAR) as a function of the number of components.
figure;
plot(1:10,cumsum(100*PCTVAR(2,:)),'-bo');
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in y');

%% compare predicted vs measured CDOM loss
horzcat(1, Q(1,:)) * Beta(:,1)     %with intercept
Q(1,:) * Beta(2:end,1)   %without intercept


%%
figure;
Xdata = 290:700;
Ydata = 230:700;
Zdata = Beta(2:end, :)';
contourf(Xdata, Ydata, Zdata, 100, 'linecolor', 'none')
colorbar()
xlabel('Irradiance Wavelength (nm)')
ylabel('Response Wavelength (nm)')

figure;
plot([290:700], Beta(2:end, :))
%%
save('C:\NTU\Research\Photodegradation experiment\manuscript\Frontiers in Marine Science\Review and Revisions\PLSR.mat');



%% now do gradient descent
load('C:\NTU\Research\Photodegradation experiment\manuscript\Frontiers in Marine Science\Review and Revisions\PLSR.mat');


%% do it for a single response wavelength each time
Beta_opt = zeros(411,471) .* nan;
J_history_all = zeros(2000,471) .* nan;

X = Q;
alpha = 100;
num_iters = 2000;


for i = 1:471
    y = CDOMloss(:,i);
    theta = Beta(2:end,i);

    [theta_opt, J_history] = gradientDescent(X, y, theta, alpha, num_iters);
    Beta_opt(:,i) = theta_opt;
    J_history_all(:,i) = J_history;
end


%%
figure;
plot(290:700, Beta_opt(:,1),'r-');hold on
plot(290:700, Beta(2:end,1),'k-')

%%
figure('color','w');
Xdata = [290:700];
Ydata = [230:700];
Zdata = Beta_opt';
contourf(Xdata, Ydata, Zdata, 100, 'linecolor', 'none')
a = colorbar();
a.Label.String = 'AQY (m^2 mol-photons^-^1 nm^-^1';
xlabel('Irradiance Wavelength (nm)')
ylabel('Response Wavelength (nm)')
set(gca, 'fontsize', 15)


%% Do some comparisons

wl = [230:700];

Pred = X * Beta_opt;
for k = 1:11
    plot(wl, CDOMloss(k,:), 'k-');hold on
    plot(wl, Pred(k,:), 'k--');hold off
    pause(1)
    
end


%% make compiled plots for comparison
wl = [230:700];
Pred = X * Beta_opt;

figure('color','w');
for k = 1:11
subplot(3,4,k)
plot(wl, CDOMloss(k,:), 'k-');hold on
plot(wl, Pred(k,:), 'k--');
legend('Measured','Predicted')
title(['Cut-off ',num2str(Cut_off(k)), ' nm'])
xlabel('CDOM Wavelength (nm)')
ylabel('△a_g (m^2)')
end


%%
save('C:\NTU\Research\Photodegradation experiment\manuscript\Frontiers in Marine Science\Review and Revisions\CDOM_AQY_Zhu.mat')





%% functions
function [theta_opt, J_history] = gradientDescent(X, y, theta, alpha, num_iters)
J_history = zeros(num_iters, 1); %record the decreasing of cost function

for iter = 1:num_iters % upper limit for the number of iteration
    [J_history(iter), grad]= computeCost(X, y, theta); % compute the value of cost
    %function and its partial derivatives
    theta = theta - alpha * grad; % Equation S3-6
end

theta_opt = theta;
end


%%
function [J,grad] = computeCost(X, y, theta)
  m = size(y,1);
  y_pred = sum(X .* theta', 2);   %calculate the predicted CDOM absorption loss using the theta
  J = 0.5/m .* sum((y_pred - y).^2, 1);   %cost function, mean squared error between predicted and measured loss of absorption coefficient
  grad = 1/m .* sum((y_pred - y) .* X, 1);   %gradient along the theta vector
  grad = grad(:);   %turn the matrix to a single column
end
% 