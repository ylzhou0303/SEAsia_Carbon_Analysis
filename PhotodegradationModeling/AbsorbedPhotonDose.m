% import data
month = 7;

% pick the appropriate solar spectra according to which month you are
% in

E = E_all(:,:,month)';
SZA = SZA_all(:,:,month)';

% calculate the Kd
CDOMforKd = CDOM_t(1, 51:end);  % update the Kd using the CDOM at each day
A = CDOMforKd + meanPABS.pabs' + W.abs';  % the total absorption is the sum of CDOM, PABS and water abs
Kd = (1 + 0.005 .* SZA) .* A + 4.18 .* (1 - 0.52 .* exp(-10.8 * A)) .* BB; %Lee et al.2005

%% calculate the absorbed photons within 0 - 0.2m
z = 0.05;
thickness = 0.1;

Q_surface = E .* S .* exp(-Kd .* z) .* CDOMforKd .* thickness;


%% calculate the absorbed photons for the entire water column
F_CDOM = ((1 + 0.005 .* SZA) .* CDOM_t(1,51:end)) ./ Kd; %only take 300-700nm of the CDOM_t to match the length of Kd

Q_total = E .* S .* (1 - exp(-Kd .* MLD)) .* F_CDOM;



%%
figure('color','w')
wl = 300:700;
plot(wl, Q_surface(6,:), 'b-');hold on
plot(wl, Q_total(6,:), 'k-')

legend('0-0.1m','Entire water column')
ylabel('Absorbed photons (mol s^-^1 nm^-^1)')
xlabel('Wavelength(nm)')


%%
Q_integrated = sum(Q_total(6,:)); %the total number of photons absorbed across the water column

Q_joey = Q_surface(6,:) .* 7.6137;   % this spectrum has the same shape as the surface, but same integrated absorbed photon dose as the entire wate column

%%
figure('color','w')
plot(wl, Q_total(6,:), 'k-');hold on
plot(wl, Q_surface(6,:),'b-')
plot(wl, Q_joey, 'b--');

legend('Entire water column','Proportionally increased suface','Surface, 0-0.1m')
ylabel('Absorbed photons (mol s^-^1 nm^-^1)')
xlabel('Wavelength(nm)')

%% import AQY from Exp 4
AQY_spc = 0.0022 * exp(-0.0086 .* wl);
DOC_loss1 = sum(AQY_spc .* Q_total(6,:));
DOC_loss2 = sum(AQY_spc .* Q_joey);

DOC_loss3 = 69 .* 10^-6 .* Q_integrated;
