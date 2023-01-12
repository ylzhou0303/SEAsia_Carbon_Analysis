% COMPILATION OF ALL THE PARAMETERS NEEDED FOR THE MODELLING: PABS, bb,
% AQY, E(scalar, below surface), DOC, CDOM
% wavelength range: 300 - 700nm (PABS measurement starts from 300nm)
% It means that the light between 290 and 300nm will have to be neglected for the
% modelling


%% import the AQY for CDOM
clear
load('F:\NTU\Research\Photodegradation experiment\photodeg modelling\MeanCDOMAQY.mat','MeanCDOM_AQY','MeanCDOM_AQY_std')

AQY_CDOM_250_700 = table();
AQY_CDOM_250_700.wl = [250:700]';   % model the change in CDOM spectrum of 250-700nm. I can only start from 250nm, because the CDOM data of Sematan, Samunsam starts from 250nm, due to the NaN3.
SolarLight_wl = [300:700];   % For the AQY spectrum, it should match the solar light spectrum
AQY_CDOM_250_700{:,2:402} = MeanCDOM_AQY(21:end,:);


%% import AQY for DOC loss
load('F:\NTU\Research\Photodegradation experiment\photodeg modelling\AQYs_DOC.mat');

%% for Talang 
load('F:\NTU\Research\Photodegradation experiment\Talang, Sarawak PABS bb\Talang_PABSfit_bb.mat')
meanPABS = Talang_MeanPabs; % I have already calculated the mean and std for Talang PABS and bb
meanbb = Talang_Meanbb;

load('F:\NTU\Research\Photodegradation experiment\Kd\Kdmodelling.mat','W')
W = W(101:501,:);

%% E (below surface, scalar)
location = 'Talang';
condition = 'clear';  % need to change the condition here if needed

matname = ['Spectra_od_',location,'_',condition,'.mat'];
load(['F:\NTU\Research\Photodegradation experiment\TUV model\TUV outputs\',matname],'Spectra_od','SZA');

SolarSpc = table();
SolarSpc.UTC = SZA.UTC;
SolarSpc.SZA = SZA.sza;

Spectra_tailored = Spectra_od(:,71:471);

% convert to Einst
plank = 6.62607004*10^-34; %Plank constant
c = 299792458;  % speed of light   m/s
avog = 6.022*10^23; % Avogadro constant
wl = [300:700];

Spc_NumPho = Spectra_tailored .* ( wl .* 10^-9) ./(plank * c);  %convert energy to number of photons   n m-2 s-1
Spc_Einst = Spc_NumPho ./ avog;       % convert from number of photons to mol of photons   mol m-2 s-1
SolarSpc{:,3:403} = Spc_Einst;  % scalar solar spectrum below surface, unit: mol photon m-2 nm-1 s-1, range: 300-700 nm


%% DOC and CDOM at T0 for the Talang region
% The riverine DOC and CDOM takes the mean of Samunsam and Sematan River
%DOC_river = mean([239.6 1188.0 1798.7 402.8]); %unit: umol L-1
DOC_river = mean([1188.0 1798.7]); %unit: umol L-1
DOC_t0 = DOC_river .* (1 - 29/33);


%CDOM_Sematan1 = [14.9194	14.9786	14.991	14.9918	14.962	14.896	14.8242	14.7496	14.67	14.5674	14.4742	14.4114	14.3342	14.211	14.1048	13.751	13.6314	13.521	13.4162	13.3018	13.1924	13.0534	12.9254	12.8218	12.7064	12.587	12.4544	12.3314	12.2142	12.0936	11.9636	11.8324	11.7014	11.5822	11.4612	11.3234	11.199	11.0818	10.9758	10.8452	10.7136	10.5908	10.4806	10.3622	10.2306	10.111	9.9984	9.875	9.7724	9.6778	9.5686	9.4604	9.3494	9.2492	9.1432	9.0446	8.9402	8.8378	8.7324	8.627	8.5318	8.4342	8.3404	8.2494	8.147	8.0458	7.964	7.8714	7.7744	7.679	7.5926	7.5088	7.4144	7.3198	7.229	7.1476	7.064	6.9766	6.893	6.8044	6.7194	6.6342	6.5486	6.4636	6.3844	6.3044	6.2208	6.1446	6.0634	5.988	5.9098	5.8292	5.7498	5.6748	5.5952	5.521	5.4514	5.3796	5.3032	5.227	5.1548	5.0834	5.0166	4.9436	4.8724	4.805	4.7366	4.668	4.5992	4.5424	4.4702	4.3998	4.3334	4.2812	4.222	4.1526	4.09	4.031	3.9718	3.9118	3.8502	3.79	3.7296	3.6736	3.616	3.5576	3.503	3.4474	3.3974	3.3492	3.2924	3.24	3.1914	3.14	3.089	3.0418	2.994	2.9458	2.9008	2.8562	2.8096	2.7626	2.72	2.681	2.637	2.5986	2.561	2.519	2.4788	2.4422	2.4038	2.3666	2.3296	2.2932	2.2568	2.2226	2.191	2.1558	2.1246	2.0894	2.0574	2.027	2.0004	1.9658	1.9368	1.906	1.8772	1.8488	1.8198	1.7924	1.7654	1.7388	1.716	1.691	1.6644	1.6394	1.6176	1.593	1.5718	1.5484	1.5258	1.5014	1.48	1.4584	1.4382	1.4182	1.3978	1.3752	1.3566	1.3364	1.317	1.299	1.2804	1.2648	1.2464	1.2278	1.2122	1.195	1.1814	1.1654	1.15	1.1328	1.1174	1.1014	1.0896	1.074	1.061	1.0474	1.0334	1.0198	1.0068	0.9958	0.986	0.9724	0.961	0.949	0.9368	0.9256	0.9156	0.9048	0.8934	0.8826	0.8724	0.8626	0.8532	0.843	0.8332	0.8228	0.8138	0.8044	0.7964	0.787	0.7772	0.7668	0.7582	0.752	0.745	0.7346	0.7276	0.721	0.7132	0.7036	0.6976	0.6924	0.6844	0.6772	0.6684	0.6604	0.6548	0.6472	0.6406	0.6328	0.6272	0.6194	0.6122	0.6054	0.6002	0.5936	0.5872	0.5812	0.5752	0.5682	0.5628	0.5568	0.5498	0.544	0.539	0.5324	0.5258	0.5214	0.5162	0.51	0.5046	0.4996	0.4944	0.4886	0.484	0.477	0.4722	0.4672	0.4632	0.4594	0.4542	0.4476	0.4448	0.4396	0.4362	0.4302	0.4258	0.4234	0.4182	0.4134	0.4088	0.4058	0.4028	0.3972	0.392	0.3878	0.385	0.3802	0.3762	0.3726	0.3686	0.3654	0.362	0.3582	0.3538	0.3512	0.3482	0.3444	0.3422	0.3386	0.3346	0.331	0.328	0.3234	0.3202	0.3192	0.315	0.3112	0.309	0.306	0.3018	0.3006	0.2966	0.2942	0.2914	0.287	0.2858	0.2832	0.2806	0.2784	0.274	0.273	0.2686	0.2678	0.2638	0.262	0.2594	0.2564	0.254	0.2518	0.25	0.2486	0.2454	0.2428	0.24	0.2384	0.2354	0.2344	0.2304	0.2298	0.2264	0.2248	0.2214	0.2202	0.2196	0.2186	0.2148	0.213	0.2122	0.209	0.208	0.206	0.2048	0.2036	0.1996	0.1994	0.1938	0.1928	0.1924	0.1888	0.1874	0.1872	0.1842	0.1838	0.1814	0.18	0.1792	0.1766	0.1756	0.1728	0.1718	0.1698	0.1682	0.1662	0.1656	0.1634	0.1636	0.1612	0.1578	0.159	0.1556	0.156	0.1542	0.1512	0.1506	0.1496	0.1474	0.146	0.1462	0.1444	0.142	0.142	0.14	0.1354	0.1388	0.135	0.1344	0.132	0.1308	0.1304	0.13	0.1284	0.1284	0.1248	0.1292	0.1238	0.1232	0.1216	0.12	0.1196	0.1202	0.1158	0.117	0.1136	0.114	0.1118	0.1124	0.1132	0.1092	0.108	0.1062	0.1064	0.1058	0.104	0.1056	0.1008	0.1022	0.1014	0.0982	0.0998	0.0986	0.0956	0.0954	0.0944	0.0924	0.094	0.0908	0.0918	0.0902	0.0886	0.0878	0.0868	0.0866	0.0854	0.0862	0.0846	0.0846	0.0852	0.0824	0.0828	0.0818	0.0796	0.0816	0.0788	0.079	0.0786	0.0784	0.0756	0.0766	0.0772	0.0756	0.0736	0.0746	0.0722	0.0708	0.0698	0.0714	0.069	0.0702	0.0676	0.0676	0.067	0.0662	0.065	0.0654	0.0638	0.0638	0.063	0.064	0.06	0.0618	0.0582	0.0616	0.0586	0.06	0.0568	0.0584	0.0608	0.0584	0.0586	0.0584	0.0548	0.0572	0.0568	0.0568	0.0548	0.0552	0.0568	0.0574	0.0542	0.0568	0.053	0.0534	0.0542	0.0536	0.0546	0.0524	0.055	0.0516	0.05	0.0548	0.0528	0.0498	0.049	0.0492	0.051	0.0522	0.0492	0.049	0.0464	0.0494	0.0492	0.0498	0.0496	0.0452	0.0456	0.0468	0.048	0.0468	0.0468	0.0486	0.0458	0.0454	0.0454	0.0436	0.0474	0.0422	0.0474	0.043	0.045	0.043	0.0412	0.0442	0.0442	0.044	0.0428	0.044	0.0452	0.0406	0.0428	0.0398	0.0422	0.0394	0.037	0.0402	0.0394	0.0386	0.039	0.039	0.0376	0.041	0.0408	0.039	0.0388	0.0392	0.039	0.0366	0.0374	0.038	0.036	0.0338	0.0364	0.035	0.035	0.0368	0.0362	0.0348	0.0348	0.037	0.0328	0.0364	0.0388	0.0322	0.0334	0.0296	0.0282	0.0316	0.032	0.0308	0.037	0.0356	0.0316	0.0294	0.0364	0.032	0.0294	0.0306	0.0374	0.0368	0.0334	0.0332	0.0326	0.0324	0.033	0.0352	0.0314	0.0334	0.0336	0.0284	0.0344	0.0316	0.033	0.0328	0.0308	0.035	0.0292	0.0324	0.0328	0.0334	0.0318	0.0322	0.033	0.0328	0.034	0.0308	0.0318	0.036	0.0306	0.032	0.0306	0.0336	0.0322	0.034	0.0322	0.0302	0.032	0.03	0.0322	0.026];
CDOM_Samunsam1 = [84.056	83.632	83.242	82.808	82.372	81.93	81.462	81.014	80.502	79.912	79.38	78.92	78.552	78.082	77.63	76.54	75.874	75.26	74.684	74.098	73.512	72.73	72.02	71.47	70.942	70.316	69.644	68.94	68.374	67.736	67.066	66.342	65.63	64.998	64.37	63.664	62.984	62.378	61.83	61.2	60.508	59.856	59.268	58.722	57.992	57.386	56.746	56.082	55.536	55.06	54.536	53.954	53.4	52.814	52.288	51.752	51.172	50.586	50.012	49.42	48.914	48.366	47.876	47.364	46.82	46.252	45.852	45.314	44.806	44.276	43.818	43.336	42.844	42.334	41.87	41.446	41.014	40.562	40.076	39.644	39.192	38.736	38.3	37.846	37.424	37.004	36.55	36.186	35.736	35.302	34.908	34.486	34.088	33.646	33.228	32.802	32.43	32.036	31.64	31.186	30.838	30.42	30.06	29.618	29.25	28.872	28.498	28.096	27.734	27.376	26.996	26.582	26.236	25.904	25.624	25.194	24.846	24.504	24.16	23.844	23.52	23.158	22.794	22.506	22.15	21.842	21.496	21.198	20.906	20.62	20.322	20.026	19.702	19.426	19.096	18.842	18.51	18.266	17.994	17.732	17.438	17.222	16.932	16.716	16.456	16.224	16.004	15.734	15.516	15.256	15.038	14.85	14.608	14.374	14.144	13.94	13.764	13.57	13.37	13.18	12.962	12.782	12.646	12.482	12.246	12.078	11.906	11.732	11.544	11.402	11.258	11.09	10.988	10.83	10.658	10.518	10.384	10.244	10.114	9.992	9.848	9.718	9.564	9.46	9.338	9.208	9.074	8.956	8.858	8.75	8.61	8.544	8.438	8.348	8.228	8.146	8.012	7.938	7.87	7.762	7.654	7.572	7.486	7.396	7.326	7.226	7.168	7.07	6.974	6.926	6.808	6.754	6.692	6.612	6.534	6.504	6.406	6.324	6.324	6.21	6.148	6.086	6.022	5.97	5.882	5.854	5.768	5.72	5.65	5.594	5.532	5.504	5.412	5.33	5.298	5.274	5.222	5.182	5.11	5.086	5.052	4.97	4.928	4.892	4.87	4.812	4.758	4.748	4.662	4.638	4.608	4.554	4.504	4.458	4.428	4.364	4.36	4.266	4.264	4.216	4.194	4.152	4.104	4.082	4.032	3.978	3.948	3.902	3.89	3.85	3.818	3.762	3.716	3.706	3.674	3.598	3.576	3.554	3.506	3.444	3.438	3.4	3.352	3.322	3.298	3.282	3.266	3.218	3.178	3.146	3.118	3.088	3.022	3.026	3.01	2.956	2.946	2.912	2.88	2.864	2.842	2.808	2.778	2.75	2.726	2.688	2.652	2.64	2.626	2.594	2.582	2.536	2.534	2.5	2.482	2.462	2.418	2.384	2.39	2.344	2.358	2.308	2.282	2.296	2.234	2.23	2.212	2.184	2.154	2.13	2.114	2.086	2.066	2.056	2.026	2.016	1.99	1.982	1.95	1.924	1.914	1.914	1.894	1.866	1.838	1.852	1.788	1.786	1.784	1.788	1.754	1.742	1.706	1.672	1.672	1.67	1.644	1.63	1.602	1.592	1.584	1.574	1.552	1.554	1.512	1.522	1.474	1.468	1.45	1.436	1.404	1.41	1.408	1.396	1.386	1.356	1.352	1.33	1.328	1.322	1.284	1.274	1.25	1.244	1.246	1.222	1.216	1.208	1.178	1.16	1.146	1.144	1.126	1.106	1.128	1.114	1.094	1.062	1.068	1.064	1.062	1.04	1.014	0.986	1.016	0.988	0.96	0.964	0.938	0.94	0.932	0.922	0.924	0.896	0.888	0.884	0.87	0.858	0.826	0.83	0.856	0.816	0.802	0.796	0.812	0.762	0.796	0.766	0.762	0.758	0.74	0.734	0.732	0.73	0.73	0.722	0.714	0.668	0.696	0.67	0.666	0.65	0.654	0.66	0.656	0.624	0.62	0.628	0.594	0.6	0.592	0.602	0.598	0.602	0.578	0.554	0.538	0.552	0.55	0.554	0.536	0.526	0.54	0.51	0.504	0.49	0.498	0.5	0.478	0.452	0.49	0.476	0.446	0.472	0.442	0.462	0.43	0.456	0.43	0.446	0.42	0.416	0.422	0.388	0.4	0.368	0.406	0.382	0.38	0.37	0.378	0.36	0.374	0.382	0.398	0.34	0.348	0.354	0.352	0.33	0.362	0.328	0.314	0.318	0.342	0.352	0.312	0.3	0.328	0.312	0.286	0.312	0.324	0.286	0.284	0.29	0.28	0.286	0.278	0.294	0.284	0.278	0.322	0.256	0.266	0.242	0.268	0.29	0.248	0.238	0.232	0.246	0.224	0.236	0.232	0.218	0.27	0.24	0.24	0.238	0.23	0.23	0.236	0.206	0.196	0.246	0.21	0.202	0.206	0.204	0.202	0.218	0.186	0.192	0.162	0.19	0.192	0.168	0.188	0.16	0.2	0.168	0.178	0.196	0.202	0.17	0.18	0.13	0.134	0.162	0.166	0.174	0.138	0.15	0.16	0.136	0.154	0.168	0.138	0.142	0.128	0.106	0.126	0.134	0.12	0.102	0.112	0.114	0.118	0.11	0.12	0.116	0.106	0.102	0.064	0.12	0.086	0.088	0.114	0.09	0.106	0.062	0.088	0.096	0.07	0.044	0.086	0.086	0.12	0.072	0.118	0.096	0.056	0.066	0.092	0.09	0.074	0.042	0.07	0.116	0.048	0.11	0.072	0.094	0.04	0.082	0.076	0.05	0.13	0.11	0.112	0.066	0.082	0.066	0.08	0.114	0.072	0.09	0.122	0.13	0.108	0.088	0.112	0.098	0.07	0.08	0.108	0.066	0.068	0.088	0.114];
%CDOM_Sematan2 = [20.5000	20.3900	20.1600	20.0300	19.9500	19.8100	19.6400	19.4400	19.2400	18.8505	18.5210	18.3315	18.1820	17.9625	17.7530	17.3835	17.1340	16.8945	16.7350	16.5655	16.3860	16.1665	15.9770	15.7975	15.6480	15.5060	15.3420	15.1700	15.0060	14.8520	14.6760	14.4940	14.2860	14.1380	13.9760	13.7680	13.5980	13.4280	13.2760	13.1240	12.9440	12.7600	12.5960	12.4460	12.2560	12.1080	11.9500	11.7820	11.6000	11.5060	11.3540	11.2200	11.0760	10.9780	10.8220	10.6900	10.5600	10.4040	10.2820	10.1400	10.0260	9.9160	9.7840	9.6560	9.5700	9.4180	9.2880	9.1800	9.0480	8.9360	8.8260	8.7180	8.5880	8.4800	8.3540	8.2500	8.1280	8.0120	7.9220	7.8100	7.6820	7.5720	7.4580	7.3520	7.2660	7.1420	7.0440	6.9500	6.8500	6.7320	6.6340	6.5220	6.4580	6.3220	6.2360	6.1660	6.0680	5.9580	5.8700	5.7520	5.6720	5.6220	5.5200	5.4380	5.3380	5.2560	5.1780	5.0960	5.0080	4.9440	4.8900	4.7880	4.6780	4.6520	4.5800	4.5060	4.4400	4.3920	4.2860	4.2240	4.1780	4.1120	4.0140	3.9920	3.8940	3.8360	3.7800	3.7000	3.6740	3.5800	3.5220	3.4760	3.4220	3.3180	3.2760	3.2400	3.2120	3.1360	3.0820	3.0140	2.9880	2.9540	2.8980	2.8460	2.8040	2.7620	2.7100	2.6760	2.6180	2.5600	2.5340	2.5040	2.4280	2.3960	2.3660	2.3420	2.3040	2.2560	2.2420	2.1840	2.1560	2.1160	2.0960	2.0700	2.0360	1.9720	1.9560	1.9520	1.9100	1.8760	1.8340	1.8180	1.7940	1.7320	1.7380	1.7040	1.6500	1.6640	1.6400	1.6340	1.6000	1.5500	1.5280	1.5380	1.4860	1.4660	1.4400	1.4380	1.4460	1.4040	1.3700	1.3580	1.3380	1.3560	1.3200	1.2900	1.2800	1.2940	1.2600	1.2040	1.2120	1.1660	1.1820	1.1860	1.1540	1.1380	1.1360	1.0880	1.1020	1.1000	1.0580	1.0420	1.0360	1.0320	1.0140	1.0060	1.0080	0.9820	0.9920	0.9940	0.9400	0.9640	0.9200	0.9220	0.8900	0.8840	0.8680	0.8620	0.8460	0.8360	0.8340	0.8120	0.7880	0.7980	0.7720	0.7700	0.7600	0.7180	0.7320	0.7180	0.7120	0.7000	0.6880	0.7180	0.6900	0.6560	0.7020	0.6600	0.6520	0.6420	0.6300	0.6320	0.6380	0.6140	0.5980	0.6140	0.6020	0.6020	0.5720	0.5660	0.5960	0.5600	0.5640	0.5500	0.5660	0.5360	0.5260	0.5180	0.5100	0.5560	0.5180	0.4920	0.5020	0.4880	0.5020	0.4900	0.4920	0.4720	0.4940	0.4360	0.4800	0.4580	0.4380	0.4460	0.4600	0.4100	0.4340	0.4460	0.4520	0.4020	0.4140	0.4240	0.4100	0.3940	0.3940	0.3860	0.3980	0.3700	0.3900	0.3800	0.3860	0.3420	0.3820	0.3620	0.3600	0.3600	0.3540	0.3560	0.3560	0.3540	0.3480	0.3300	0.3240	0.3340	0.3400	0.3220	0.2980	0.3280	0.2820	0.3280	0.3080	0.2820	0.2860	0.3060	0.3020	0.3140	0.2840	0.2860	0.3060	0.2940	0.2980	0.2620	0.3020	0.2620	0.2860	0.2760	0.2720	0.2560	0.2660	0.2800	0.2720	0.3020	0.2500	0.2620	0.2740	0.2500	0.2520	0.2080	0.2320	0.2200	0.2560	0.2560	0.2020	0.2200	0.2500	0.2140	0.2600	0.2400	0.1920	0.2080	0.2160	0.2140	0.2160	0.2300	0.2020	0.1800	0.1980	0.1880	0.1820	0.2000	0.1940	0.1940	0.1780	0.1920	0.2060	0.1700	0.1740	0.1900	0.1860	0.1800	0.1760	0.1940	0.1680	0.1660	0.1920	0.1920	0.1260	0.1880	0.1520	0.1480	0.1660	0.2240	0.1980	0.1500	0.1580	0.1700	0.1600	0.1420	0.1500	0.1560	0.1620	0.1680	0.1540	0.1300	0.1440	0.1700	0.1520	0.1460	0.1240	0.1320	0.1300	0.1180	0.1440	0.1300	0.1420	0.1080	0.1360	0.1400	0.1380	0.1480	0.1460	0.1120	0.1120	0.1020	0.1420	0.1180	0.1120	0.1160	0.0920	0.1080	0.1220	0.1420	0.1240	0.1140	0.1080	0.0860	0.1080	0.1080	0.0940	0.0900	0.1440	0.1020	0.0760	0.0780	0.0860	0.1020	0.0920	0.1140	0.0780	0.0640	0.0880	0.0540	0.0980	0.0820	0.0560	0.0920	0.0880	0.0820	0.0740	0.0820	0.0620	0.0560	0.0820	0.1040	0.0720	0.0720	0.0760	0.0480	0.0780	0.0820	0.0600	0.0680	0.0720	0.0760	0.0720	0.0760	0.0540	0.0520	0.0760	0.0680	0.0760	0.0700	0.0740	0.0280	0.0900	0.0200	0.0640	0.0620	0.0700	0.0580	0.0780	0.0720	0.0540	0.0580	0.0620	0.0500	0.0480	0.0780	0.0460	0.0660	0.0740	0.0600	0.0560	0.0480	0.0720	0.0620	0.0720	0.0520	0.0540	0.0900	0.0600	0.0700	0.0240	0.0420	0.0580	0.0780	0.0580	0.0680	0.0580	0.0760	0.0840	0.0740	0.0660	0.0540	0.0700	0.0480	0.0160	0.0440	0.0640	0.0600	0.0560	0.1100	0.0700	0.0660	0.0360	0.0900	0.0520	0.0660	0.0420	0.0760	0.0340	0.0520	0.0440	0.0840	0.0700	0.1160	0.0140	0.0320	0.0720	0.0520	0.0600	0.0160	0.0980	0.0660	0.1080	0.0860	0.0880	0.0320	0.0620	0.0660	0.0400	0.0580	0.0360	0.0420	0.0060	0.0920	0.0540	0.0100	0.0180	0.0340	0.0560	0.0500	0.0400	0.0380	0.0640	0.0640	0.0340	0.0620	0.0600	0.0480	0.0500	0.0600	0.0500	0.0660	0.0740	0.0440	0.0520	0.0260	0.0520	0.0300	0.0120	0.0280	-0.0100	0.0360	-0.0040	-0.0220	0.0280	0.0060	0.0580	0.0620	0.0580	0.0160	0.0320	0.0440	0.0260	0.0440	0.0500	0.0720	0.0380	0.0720	0.0120	0.0640	0.0360	0.0660	0.0340	0.0380	0.0440	0.0680	0.0380	0.0680	0.0180	0.0500	0.0660	0.0260	0.0540	0.0200	0.0600	0.0360	0.0540	0.0520	0.0560	0.0760	0.0520	0.0640	0.0680	0.0620	0.0880	0.0620	0.0580	0.0880	0.0360	0.0400	0.0900	0.0840	0.0820	0.0320	0.0900	0.0500	0.0760	0.0680	0.0920];
CDOM_Samunsam2 = [124.5500	124.0000	123.1500	122.4000	122.1000	121.4000	120.7500	119.8500	119.0500	117.7500	116.6000	115.8000	115.6500	115.1000	114.4500	112.6000	111.6000	110.4500	109.5875	108.6375	107.6375	106.4250	105.2500	104.4500	103.5250	102.7000	101.6125	100.5750	99.5750	98.6000	97.4875	96.3250	95.1875	94.2875	93.3500	92.2375	91.1875	90.0625	89.3250	88.4500	87.2750	86.4500	85.4125	84.5125	83.5000	82.4250	81.5375	80.5250	79.6500	78.9125	78.1375	77.2750	76.3625	75.6000	74.7625	73.8750	73.0125	72.1625	71.2625	70.3250	69.5875	68.7250	67.9250	67.2125	66.3875	65.4750	64.8000	64.1750	63.3625	62.5500	61.8750	61.2500	60.4625	59.6250	58.9250	58.3125	57.5375	56.9875	56.3500	55.6375	54.9125	54.3125	53.6875	53.0125	52.3250	51.8375	51.0875	50.4875	49.8125	49.2375	48.6750	48.0125	47.3250	46.8000	46.0375	45.4125	45.0000	44.4000	43.9000	43.1625	42.6125	41.9875	41.4375	40.8500	40.2375	39.7000	38.9875	38.4500	37.9500	37.4750	36.7875	36.2375	35.5500	35.1375	34.7250	34.1000	33.7000	33.1375	32.6000	32.0750	31.5125	31.0000	30.6125	30.1000	29.5500	29.1250	28.6375	28.2375	27.8125	27.3625	26.9000	26.5125	26.0125	25.5750	25.2125	24.8250	24.4125	24.0125	23.6250	23.2875	22.8125	22.3750	22.0875	21.7375	21.3375	21.0625	20.7875	20.3500	20.0000	19.6875	19.4000	19.1000	18.7375	18.5000	18.2000	17.9000	17.5625	17.3125	17.0125	16.7875	16.4625	16.2125	15.8625	15.7250	15.5500	15.2250	14.9500	14.6875	14.5125	14.2125	14.1375	13.8500	13.6375	13.4375	13.2375	13.0375	12.8375	12.6625	12.5125	12.3875	12.2750	12.0375	11.7625	11.5500	11.5250	11.3125	11.1000	10.9750	10.9375	10.6750	10.5250	10.5750	10.2625	10.2250	10.0750	10.0125	9.7250	9.6750	9.6125	9.4125	9.3750	9.1750	9.1250	8.8500	8.7375	8.6875	8.5750	8.5000	8.3250	8.1875	8.0375	7.9750	7.8875	7.8000	7.6875	7.6500	7.5250	7.5125	7.2750	7.2625	7.2000	7.1125	7.0250	6.9375	6.8750	6.7625	6.6875	6.5500	6.5375	6.4875	6.3875	6.2375	6.2000	6.1250	6.0625	6.0000	5.9500	5.8375	5.7750	5.8000	5.8000	5.6250	5.4750	5.6000	5.5750	5.4625	5.4250	5.2500	5.2625	5.2625	5.2500	5.1250	5.1000	5.0125	4.9375	4.9500	4.8500	4.8250	4.7625	4.6875	4.6625	4.6500	4.5625	4.5750	4.5125	4.4875	4.4375	4.2750	4.3500	4.3125	4.2125	4.2125	4.1125	4.1125	4.0375	3.9500	3.9625	3.8875	3.9500	3.7125	3.7250	3.6500	3.6375	3.4625	3.5375	3.5250	3.3750	3.4375	3.4625	3.3250	3.3250	3.1875	3.2000	3.2000	3.0625	3.1125	3.0375	3.0000	3.0125	3.0125	2.8625	2.9250	2.8875	2.8500	2.8500	2.6875	2.6250	2.7250	2.6000	2.5625	2.6250	2.6000	2.5125	2.5625	2.4125	2.4125	2.5250	2.3750	2.3625	2.4000	2.3500	2.2875	2.2000	2.1500	2.1625	2.1750	2.1000	2.1000	2.1625	2.0875	2.0875	1.9875	2.1000	1.9875	2.0125	2.0500	1.9625	1.9875	2.0625	1.8875	1.9250	1.8875	1.7750	1.8500	1.8250	1.8125	1.8375	1.7125	1.7625	1.8500	1.7250	1.6750	1.6375	1.7125	1.5500	1.6000	1.5625	1.6375	1.4750	1.5125	1.6250	1.5875	1.4125	1.5375	1.4125	1.3875	1.4000	1.4000	1.4125	1.3625	1.3125	1.3250	1.3875	1.2875	1.1750	1.2625	1.2250	1.2250	1.2250	1.2250	1.1625	1.3875	1.1750	1.0875	1.1000	1.1250	1.1250	1.0250	1.1250	1.1625	1.2250	1.1625	0.9375	1.1250	1.0250	1.0250	0.9750	0.9875	0.9125	0.9500	0.8625	0.9250	0.9125	1.0250	0.9250	0.9750	0.9500	0.9250	0.9000	0.9000	0.8625	0.8500	0.8375	0.8125	0.7750	0.8375	0.7250	0.7125	0.7500	0.7750	0.6500	0.5875	0.6500	0.8250	0.8875	0.7000	0.7500	0.5625	0.7875	0.6250	0.7375	0.6500	0.6500	0.6500	0.7000	0.6250	0.6500	0.6000	0.6375	0.6250	0.6500	0.5000	0.5500	0.6000	0.6125	0.5250	0.4750	0.5250	0.5000	0.5000	0.4250	0.4750	0.4125	0.5000	0.4625	0.3250	0.4375	0.4500	0.3750	0.4125	0.5250	0.5250	0.4125	0.3875	0.3500	0.3250	0.2000	0.3625	0.3875	0.3000	0.3250	0.2625	0.3125	0.3875	0.3750	0.3250	0.2000	0.2125	0.1750	0.3250	0.1875	0.2625	0.3250	0.2000	0.3625	0.1750	0.2250	0.0625	0.2250	0.3125	0.2250	0.2375	0.3125	0.2250	0.2000	0.1375	0.1875	0.1500	0.1250	0.1875	0.1875	0.2625	0.3250	0.1250	0.0500	0.1750	0.1625	0.1250	0.0875	0.2125	0.0875	0.2125	0.1375	0.1625	0.1625	0.1625	0.2625	0.1500	0.1375	0.0500	0.0625	0.0125	0.1125	0.1500	0.3375	0.0625	0.1875	0.1500	0.1125	0.2875	0.1000	0.1125	0.0125	0.0625	0.1875	0.1625	0.1500	0.0750	0.0000	0.1125	0.1375	0.1875	0.1500	0.0125	0.0750	0.1125	0.0000	0.1250	0.1250	0.1750	0.1625	0.0500	0.1000	0.2250	0.1625	0.0750	0.0875	0.0750	-0.0750	0.1000	0.1625	0.0375	0.0375	0.0625	0.0625	-0.0375	0.1500	0.1000	0.1500	-0.0125	0.0625	-0.0625	0.0000	-0.0750	0.0250	-0.0375	0.0250	0.0375	0.0625	-0.1500	-0.0500	-0.0250	-0.0750	0.1000	-0.0250	0.0125	0.0625	0.0375	-0.0500	-0.0875	-0.0750	-0.1125	0.0000	-0.0875	-0.1375	-0.1500	-0.0750	-0.3000	-0.2250	-0.0750	-0.2000	-0.1000	-0.0250	-0.1875	-0.0125	-0.0625	-0.1000	-0.1000	-0.2375	-0.3000	-0.1500	0.0375	0.0125	-0.2000	-0.3000	-0.2250	-0.1750	-0.2250	-0.0750	-0.0875	-0.1625	-0.0875	-0.0625	-0.0875	-0.0875	-0.0500	-0.1250	-0.1375	-0.0125	0.0875	-0.0500	-0.0250	-0.1500	-0.1625	-0.2000	0.0000	-0.2250	-0.0625	0.0500	-0.0750	-0.1500	0.0125	0.0500	-0.0375	-0.0125	-0.1250	-0.1500	-0.1375	0.0625	-0.1125	0.0375];

CDOM_river = mean(vertcat(CDOM_Samunsam1,CDOM_Samunsam2),1);
CDOM_t0 = table();
CDOM_t0.wl = [250:700]';
CDOM_t0.river = CDOM_river(1:451)';
CDOM_t0.shelf = CDOM_t0.river .* 2.303 .* (1 - 29/33);


%% Talang
MLD = 7.6;
MLD_std = 4.5;



%%
savename = ['PARs_',location,'_',condition,'.mat'];
clearvars -except AQY_CDOM_250_700 AQYs_DOC SolarSpc DOC_t0 CDOM_t0 meanPABS meanbb W MLD MLD_std DOC_river MeanCDOM_AQY MeanCDOM_AQY_std savename
save(['F:\NTU\Research\Photodegradation experiment\photodeg modelling\',savename]);