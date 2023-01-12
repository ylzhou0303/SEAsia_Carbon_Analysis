load('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\KusuHantuData.mat');


%% calculate marine endmembers in 2019
dic_mar = mean(vertcat(Kusu.DIC(33:38),Hantu.DIC(31:34))); % use data between 25 Feb 2019 to 15 Mar 2019
ta_mar = mean(vertcat(Kusu.TA(33:38),Hantu.TA(31:34)));
d13cdic_mar = mean(vertcat(Kusu.d13cDIC(33:38),Hantu.d13cDIC(31:34)));
sal_mar = mean(vertcat(Kusu.salinity(33:38),Hantu.salinity(31:34)));
d13cDOC_mar = mean(vertcat(Kusu.d13cDOC(33:38),Hantu.d13cDOC(31:34)));

%% the marine endmember in 2020
dic_mar = mean(vertcat(Kusu.DIC(56:57),Hantu.DIC(52)));
ta_mar = mean(vertcat(Kusu.TA(56:57),Hantu.TA(52)));
d13cdic_mar = mean(vertcat(Kusu.d13cDIC(56:57),Hantu.d13cDIC(52)),'omitnan');
sal_mar = mean(vertcat(Kusu.salinity(56:57),Hantu.salinity(52)));
d13cDOC_mar = mean(vertcat(Kusu.d13cDOC(56:57),Hantu.d13cDOC(52)));

%% the marine endmember using mean of 2019 and 2020 data
dic_mar = mean(vertcat(Kusu.DIC([33:36 56 57]),Hantu.DIC([31:34 52]))); % use data between 25 Feb 2019 to 15 Mar 2019, and late Feb and March in 2020
ta_mar = mean(vertcat(Kusu.TA([33:36 56 57]),Hantu.TA([31:34 52])));
d13cdic_mar = mean(vertcat(Kusu.d13cDIC([33:36 56 57]),Hantu.d13cDIC([31:34])),'omitnan');
sal_mar = mean(vertcat(Kusu.salinity([33:36 56 57]),Hantu.salinity([31:34 52])));
d13cDOC_mar = mean(vertcat(Kusu.d13cDOC([33:36 56 57]),Hantu.d13cDOC([31:34 52])));

%% standard deviation
dic_mar_std = std(vertcat(Kusu.DIC([33:36 56 57]),Hantu.DIC([31:34 52]))); % use data between 25 Feb 2019 to 15 Mar 2019, and late Feb and March in 2020
ta_mar_std = std(vertcat(Kusu.TA([33:36 56 57]),Hantu.TA([31:34 52])));
d13cdic_mar_std = std(vertcat(Kusu.d13cDIC([33:36 56 57]),Hantu.d13cDIC([31:34])),'omitnan');
sal_mar_std = std(vertcat(Kusu.salinity([33:36 56 57]),Hantu.salinity([31:34 52])));
d13cDOC_mar_std = std(vertcat(Kusu.d13cDOC([33:36 56 57]),Hantu.d13cDOC([31:34 52])));

%%
model = struct('dic_riv',453, 'ta_riv', 310,  'dic_mar',dic_mar, 'ta_mar', ta_mar, 'sal_mar',sal_mar, 'sal', [0:0.1:35], 'R_vpdb',  0.01123720, 'd13c_riv',-15.32, 'd13c_peat', -29, 'd13c_mar' , d13cdic_mar, 'd13cDOC_mar', d13cDOC_mar, 'd13c_tDIC', -32, 'd13c_carb', 0);
model_std = struct('d13cdic_mar_std',d13cdic_mar_std, 'd13cDOC_mar_std',d13cDOC_mar_std, 'dic_mar_std',dic_mar_std, 'sal_mar_std', sal_mar_std, 'ta_mar_std', ta_mar_std);
save('F:\NTU\Research\Kusu Hantu Biogeochem\manuscript\JGR\figures\endmember.mat','model','model_std');
