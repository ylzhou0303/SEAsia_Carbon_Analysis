% load 1cm corrected dataset
cd('F:\NTU\Research\Sarawak FDOM analysis_new\1cm\Corrected EEMs')
load CorrectedDataset_1cm.mat

% load 1cm diluted corrected dataset
cd('F:\NTU\Research\Sarawak FDOM analysis_new\1cm diluted\Corrected EEMs')
load CorrectedDataset_1cm_diluted.mat

% load 3mm corrected dataset
cd('F:\NTU\Research\Sarawak FDOM analysis_new\3mm\Corrected EEMs')
load CorrectedDataset_3mm.mat

% combine the corrected dataset
%% Reconcile the sampleID for all
SampleID_All=vertcat(SampleID_1cm,SampleID_1cm_diluted,SampleID_3mm);
filelist_all = vertcat(cellstr(filelist_eem_1cm),cellstr(filelist_eem_1cm_diluted),cellstr(filelist_eem_3mm));

% Corrected EEM
XcRU2_All=vertcat(XcRU2_1cm,XcRU2_1cm_diluted,XcRU2_3mm);

%% %% Assemble dataset, smooth EEM scatter and remove outliers

%Assemble the dataset
mydata= assembledataset(XcRU2_All,Ex,Em,'RU','ID',SampleID_All,'filename',filelist_all,[]);
%eemview(mydata,'X',[3 3],23,[],'ID')
classinfo(mydata)

%% Smooth it
%Remove noisy data then smooth the EEM
SubData=subdataset(mydata,[],mydata.Em>580,mydata.Ex<250);
%eemview(SubData,'X',[3 3],1,[],'ID')
Xs = smootheem (SubData,[15 15],[15 15],[18 18],[100 15],[1 1 1 1],[],3500,0);  %set the Ram2 [100 15], 100nm above, try to remove everything at the corner


%% See what happened
eemview (Xin, 'X', [3 3],1,[],'filename',[],[],'colorbar')
%% Remove outliers
% Sematan stn20 is identified as an outlier, delete it
Xin = subdataset (Xs,{'filename','Scoast Stn 14 20170907.csv','Scoast Stn 2 20170904.csv','Scoast Stn 21 20170909.csv','Sebuyau Stn 6 20170311 df10.csv','Sematan Stn 20 20170915.csv'},[],[])

% Eliminate faulty parts of EEMs
 Xin=zap(Xin,17,[390,500],[330 340])              %function zap only remove the faulty emission scans of sample 176 (Ex 345-350)
 Xin=zap(Xin,17,[390,490],[365 375]) 

eemview(Xin,'X',[3 3],1,[],'filename',[],[],'colorbar',[])      %Here you can start from any sample you want

%% save the smoothed EEM struct

cd('F:\NTU\Research\Sarawak FDOM analysis_new\PARAFAC analysis_new')
save SmoothedEEM_All_new.mat Xin Xs

%%
cd('F:\NTU\Research\Sarawak FDOM analysis_new\PARAFAC analysis_new')
load SmoothedEEM_All_new.mat

%% Exploratory data analysis
%Xin = subdataset (Xin, [89], [], []);  %remove Scoast stn 8 from the dataset because this EEM looks weird and the redisue between the PARAFAC model and the real sample is big
Xin_2 = subdataset (Xin, [], Xin.Em>550, []);

%Search for outliers using outliertest
Test1 = outliertest (Xin_2, [2,2], 5:6, 'nonnegativity', [], 'at once') %with constraint "nonnegativity", so make sure the output is always positive

% Visulize the spectra
 comparespectra (Test1, 5:6)  % The nonnegativity constraint is necessary to obtain the reasonable spectra
% comparespectra (Test1u, 3:7)
% 
% 
% %see the correlations between components (if components are highly correlated to each other, normalise them)
 compcorrplot (Test1, 5)
 compcorrplot (Test1, 6)
% compcorrplot (Test1, 6)


%% Normalization
% Develop a dataset in which EEM intensities are normalized to unit norm,
% in order to limit the effect of dilution and highly correlated components
% which will be the dataset modelled using PARAFAC

Xpre = normeem (Xin_2)
Test1p=outliertest(Xpre,[1,1],4:7,'nonnegativity',[],'at once');
spectralloadings(Test1p,4:7)
spectralloadings(Test1,4:7)

% Compare the spectral loadings of the new model with the old model
%View correlation in the normalized set 
% compcorrplot(Test1p,3)
% compcorrplot (Test1p,4)
% compcorrplot (Test1p,5)
% compcorrplot (Test1p,6)


%%
% function loadingsandlevarages: to identify the component or wavelength
% that has unsually high influence on the models ( Wavelengths with high leverages will have more influence on spectral shapes than other wavelengths, which could distort the results.)
%loadingsandleverages(Test1p,3)
loadingsandleverages(Test1p,4)
loadingsandleverages(Test1p,5) 
loadingsandleverages(Test1p,6)  
loadingsandleverages(Test1p,7) 
loadingsandleverages(Test1p,8)   
 
%%
Xpre2 = subdataset (Xpre, [64 65  148 152],[],[] ); %For 5comp model, remove these outliers
Test1p_2 = outliertest(Xpre2,[1,1],5:7,'nonnegativity',[],'at once');

% Check the loadings and leverages again
loadingsandleverages(Test1p_2,5) %use the Xpre2 to generate the 5comp model

%Actually the phbl Maludam samples shouldn't be outliers, so I tend to
%leave them in the model, but I will remove 64 and 65, which are the
%Samunsam sample, and try to generate a model using Xpre3

Xpre3 = subdataset (Xpre, [64 65], [] ,[])
Test1p_3 = outliertest(Xpre3,[1,1],5:7,'nonnegativity',[],'at once');

loadingsandleverages(Test1p_3,6)  
loadingsandleverages(Test1p_3,7) 

%% check the residuals
eemview({Test1p_3,Test1p_3},{'X','Model5','error_residuals'},[3 3],1,[],{'ID','ID','ID'},[],[],{'colorbar','colorbar','colorbar'})  
eemview({Test1p_3,Test1p_3},{'X','Model6','error_residuals'},[3 3],1,[],{'ID','ID','ID'},[],[],{'colorbar','colorbar','colorbar'})  
eemview({Test1p_3,Test1p_3},{'X','Model7','error_residuals'},[3 3],1,[],{'ID','ID','ID'},[],[],{'colorbar','colorbar','colorbar'})  


%Compare residuals for outliertest Model 5 and 6
compare2models (Test1p, 5,6,0.05)

%We can evaluate the feasibility of the PARAFAC model by using the
%spectralloaddings function, and checking whether some components are
%actually noise
spectralloadings (Test1p_3, 5:7)

% what should be the number of components in the model? use the specsse
% function, which shows the effect of adding more component on the fit
% model, expressed as the sum of squared error (SSE)

% Smaller SSE indicates a better model. However, a small decrease in SSE indicates that the last component added might be unnesscessary or detrimental.
specsse (Test1p, 5:7)
%According to the SSE plot, 4 or 5 components should be appropriate
%% Validation Phase
%Refine the models
% Obtain least squares models 

[LSmodel5,convg5,DSit5]=randinitanal(Xpre3,5,5,'nonnegativity',1e-8);
[LSmodel6,convg6,DSit6]=randinitanal(Xpre3,6,5,'nonnegativity',1e-8);
[LSmodel7,convg7,DSit7]=randinitanal(Xpre3,7,5,'nonnegativity',1e-8);

% Visualize the spectra from the least squares model
spectralloadings (LSmodel5,5)
comparespectra (LSmodel5,5)
fingerprint (LSmodel5,5)
loadingsandleverages(LSmodel5,5)
eemview({LSmodel5,LSmodel5},{'X','Model5','error_residuals'},[2 3]);


spectralloadings (LSmodel6,6)
comparespectra (LSmodel6,6)
fingerprint (LSmodel6,6)
loadingsandleverages(LSmodel6,6)
eemview({LSmodel6,LSmodel6},{'X','Model6','error_residuals'},[2 3]);


spectralloadings (LSmodel7,7)
comparespectra (LSmodel7,7)
fingerprint (LSmodel7,7)
loadingsandleverages(LSmodel7,7)   %No. 135 has unusually high leverage(>0.25)
eemview({LSmodel7,LSmodel7},{'X','Model7','error_residuals'},[2 3]);


%% Reverse the normalisation step to recover the true scores
% The reversed model is exactly the same as the original, except that the
% scores are no longer normalised.
LSmodel5r=normeem(LSmodel5,'reverse',5)
LSmodel6r=normeem(LSmodel6,'reverse',6)
LSmodel7r=normeem(LSmodel7,'reverse',7)

%% Validation by split half analysis
% dataset will be divided into different groups
S1 = splitds (Xpre3, [], 4 ,'alternating',{[1 2],[3 4],[1 3],[2,4],[1 4],[2 3]})
% (S4C6T3: split:4, combinations: 6, tests: 3)

%Use splitanalysis to generate 5- and 6- and 7- component PARAFAC models in each
%dataset split
A1 = splitanalysis (S1, 5:6,'nonnegativity',[],[],'A1')

% preliminary validation to see which split models are the same
% splitvalidation(A1,5,[1 2;3 4;5 6]);                %not validated
% splitvalidation(A1,6,[1 2;3 4;5 6]);        %not validated
% splitvalidation(A1,7,[1 2;3 4;5 6]);        %not validated

%Improve the split models by taking the best of multiple runs
%and/or decrease the convergence criterion

%5-component model
MyRunOptions=[10 10 10 10 10 10] %specify runs for each split individually
MyCC=[1e-8 1e-8 1e-8 1e-8 1e-8 1e-8]           %CC:convergence criterion
SaveResultsAs='A1b_10-8_5comp'
A1b_5comp=splitanalysis(S1,5,'nonnegativity',MyRunOptions,MyCC,SaveResultsAs);

splitvalidation(A1b_5comp,5,[1 2;3 4;5 6]);       %this time validated

%6-component model
MyRunOptions=[5 5 20 20 5 5] %specify runs for each split individually
MyCC=[1e-8 1e-8 1e-8 1e-8 1e-8 1e-8]           %CC:convergence criterion
SaveResultsAs='A1b_10-8_6comp'
A1b_6comp=splitanalysis(S1,6,'nonnegativity',MyRunOptions,MyCC,SaveResultsAs);

splitvalidation(A1b_6comp,6,[1 2;3 4;5 6]);       %not validated 

%complete the validation for 3 and 4 component model
%val3=splitvalidation(A3,3,[1 2;3 4; 5 6],{'AB','CD','AC','BD','AD','BC'},LSmodel3r);
%val4=splitvalidation(A1,4,[1 2;3 4; 5 6],{'AB','CD','AC','BD','AD','BC'},LSmodel4r);
val5=splitvalidation(A1b_5comp,5,[1 2;3 4; 5 6],{'AB','CD','AC','BD','AD','BC'},LSmodel5r);
val6=splitvalidation(A1b_6comp,6,[1 2;3 4; 5 6],{'AB','CD','AC','BD','AD','BC'},LSmodel6r);

%Examine error residuals for Model 3 and 4
eemview({LSmodel6,LSmodel6},{'X','Model6','error_residuals'},[2 3]);
loadingsandleverages(DSit6,6,5); %Observe leverages on model in e.g. Run 5

%Compare error residuals for least squares models 5 and 6
eemview({LSmodel5,LSmodel5,LSmodel6,LSmodel6},{'X','Model5','error_residuals','X','Model6','error_residuals'},[2 3],[],[],[],[1 1 0.05 1 1 0.05],'colorbar');


spectralloadings (val5,5)
spectralloadings (val6,6)


%% View some results prior to exporting
describecomp(LSmodel6,6)

%view Fmax for a particular model
[F3,scores3]=scores2fmax(val3,3)                               % calculate the Fmax representing the maximum intensity of each component in each sample
[F4,scores4]=scores2fmax(val4,4) 

%Describe components in terms of Ex and Em maxima
summary=describecomp(val6,6);

%% Export data
% calculate the peak positions of PARAFAC components
% peak=describecomp(LSmodel6,6)

% %Basic model export, including converting to Fmax
% modelout(val6,6,'drEEM_PARAFAC_demoX.xls');
% 
% %Export metadata along with a model
% metadata4export={'cruise','site','rep','longID','ID','date','i'}
% modelout(val6,6,'drEEM_PARAFAC_demo.xls',[],metadata4export);
% 
% %Project the full dataset (Xs) on to the model in val6 to get scores for outliers
% [F6,EmSpectra,ExSpectra,Ff,P6]=modelout(val6,6,'drEEM_PARAFAC_demo.xls',Xs);

% Export


modelout(val5,5,'val5 result.xls');
[F5,B5,C5,Ff5,P5]=modelout(val5,5,'val5 result_full dataset.xlsx',Xin_2);    % expand to the full dataset, including the outliers

modelout(val6,6,'val6 result.xls');
[F6,B6,C6,Ff6,P6]=modelout(val6,6,'val6 result_full dataset.xlsx',Xin_2);

% SAVE THE WORKSPACE!!

%Look at the model fit for the projected dataset
% eemview({P6,P6},{'X','Model5_P','error_residuals'},[3 3],[],[],[],[1 1 0.05]);
% close all

