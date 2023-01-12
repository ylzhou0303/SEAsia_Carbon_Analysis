%% In this file, the spectra just below the water surface are generated at each hour

load('F:\NTU\Research\Photodegradation experiment\TUV model\TUV outputs\AllSpectra_Talang_cloud.mat');

%%
AllSpectra = AllSpectra_Talang_cloud;
SZA = SZA_Talang;
%% calculate the transmittance from the air into the ocean for the direct radiance using the Fresnel's Law
PAR = table();
PAR.UTC = SZA.UTC;
PAR.theta1 = SZA.sza;

%% 1. Use Snell's Law to calculate the refracted angle, theta2
n1 = 1;  %refractive index for air
n2 = 1.338;  %refractive index of seawater   %Austin & Halikas 1976
PAR.theta2 = asind(n1 .* sind(PAR.theta1) ./ n2);

%% 2. transmission coefficient
PAR.ts = 2 .* n1 .* cosd(PAR.theta1) ./ (n1 .* cosd(PAR.theta1) + n2 .* cosd(PAR.theta2)); % s-polarization
PAR.tp = 2 .* n1 .* cosd(PAR.theta1) ./ (n1 .* cosd(PAR.theta2) + n2 .* cosd(PAR.theta1)); % p-polarization

%% 3. Transmittance
PAR.Ts = (n2 .* cosd(PAR.theta2)) ./ (n1 .* cosd(PAR.theta1)) .* (PAR.ts .^ 2); 
PAR.Tp = (n2 .* cosd(PAR.theta2)) ./ (n1 .* cosd(PAR.theta1)) .* (PAR.tp .^ 2); 

% the natural sunlight is unpolarized light, so we should take the
% arithmatic mean of the transmittance for s- and p-polarized light
PAR.T = mean(horzcat(PAR.Ts, PAR.Tp),2);


%% The transmittance for the diffuse fraction is assumed to be 0.934 (Fichot 2010)
%% The conversion factor to scalar irradiance
PAR.scalar_dir = 1 ./ (cosd(PAR.theta2)); %direct light
PAR.scalar_diff = ones(288,1) .* 1/0.859;  %diffused light

%% calculate the total transmitted light and conversion to scalar irradiance
Spectra_od = ones(288,570) * nan;
for k = 1:288
    if PAR.theta1(k) <= 90
        %calculate the sum of the direct and diffuse irradiance
        Spectra_od(k,:) = AllSpectra(:,3,k)' .* PAR.T(k) .* PAR.scalar_dir(k) + AllSpectra(:,4,k)' .* 0.934 .* PAR.scalar_diff(k);
    else
        Spectra_od(k,:) = zeros(1,570);     %if the solar zenith angle is more than 90degrees, it means that the sun is below horizon, so no irradiance
    end
end


%%
figure;
hold on
for i = 1:24:288
    plot(Spectra_od(i,:));
    pause(1)
end

%% arrange the datetime format
for i = 1:288
    tempchar = char(SZA.UTC(i));
    DT = [tempchar(1:4),'-',tempchar(5:6),'-',tempchar(7:8),' ',tempchar(9:10),':00:00'];
    UTCTime(i,1) = datetime(DT);
end

%%
figure;
hold on
plot(AllSpectra(:,6,5+24),'k-');
plot(Spectra_od(5+24,:),'k--');
hold off

%%
Spectra_od_Talang_cloud = Spectra_od;

%%
save('F:\NTU\Research\Photodegradation experiment\TUV model\TUV outputs\Spectra_od_Talang_cloud.mat')