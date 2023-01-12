% calculate the photochemical decay constant for the UniDOM


load('C:\NTU\Research\Photodegradation experiment\photodeg modelling\model results\PhotoDegModel_SG_cloud_.mat')
monthdays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; % number of days in each month
%%
Phi = table();
Phi.DOC_start = ones(12,1).*nan;
Phi.DOC_start(1) = DOC_t0;
Phi.DOC_end(1) = DOC_t.DOC_umolL(31);
Phi.a350_start(1) = CDOM_t0.shelf(101);
Phi.a350_end(1) = CDOM_t(31,101);

for i = 2:12
    Phi.DOC_start(i) = DOC_t.DOC_umolL(sum(monthdays([1:i-1])));
    Phi.DOC_end(i) = DOC_t.DOC_umolL(sum(monthdays([1:i])));
    Phi.a350_start(i) = CDOM_t(sum(monthdays([1:i-1])),101);
    Phi.a350_end(i) = CDOM_t(sum(monthdays([1:i])),101);
end

Phi.a350_mean = mean([Phi.a350_start Phi.a350_end],2);

Phi.DOClossrate = Phi.DOC_start - Phi.DOC_end;
Phi.phi = (log(Phi.DOC_start) - log(Phi.DOC_end))./monthdays';

kuv_w = 0.12;
Phi.kuv = kuv_w + Phi.a350_mean ./ 2.303; %CONVERT THE a350 to decadic absorption coefficient!


JOEY = (1 - exp(-Phi.kuv .* MLD)) ./ Phi.kuv; % Anderson et al. 2019 Eqn. 5-7
Phi.phi_ref = Phi.phi .*MLD ./ JOEY;
%%
save('C:\NTU\Research\Photodegradation experiment\photodeg modelling\model results\SG_phi_unidom.mat','Phi')


%% Talang
load('C:\NTU\Research\Photodegradation experiment\photodeg modelling\model results\PhotoDegModel_Talang_cloud_winter.mat')

%%
Phi = table();
for i = 1:3
    Phi.DOC_start(i) = DOC_t.DOC_umolL(30*(i-1)+1);
    Phi.DOC_end(i) = DOC_t.DOC_umolL(30*i);
    Phi.a350_start(i) = CDOM_t(30*(i-1)+1,101);
    Phi.a350_end(i) = CDOM_t(30*i,101);
   
end

%%
Phi.a350_mean = mean([Phi.a350_start Phi.a350_end],2);

Phi.DOClossrate = Phi.DOC_start - Phi.DOC_end;
Phi.phi = (log(Phi.DOC_start) - log(Phi.DOC_end))./30;

kuv_w = 0.12;
Phi.kuv = kuv_w + Phi.a350_mean ./ 2.303; %CONVERT THE a350 to decadic absorption coefficient!


JOEY = (1 - exp(-Phi.kuv .* MLD)) ./ Phi.kuv; % Anderson et al. 2019 Eqn. 5-7
Phi.phi_ref = Phi.phi .*MLD ./ JOEY;


%%
save('C:\NTU\Research\Photodegradation experiment\photodeg modelling\model results\Talang_phi_unidom.mat','Phi')