demopath = 'F:\NTU\Research\Sarawak EEM\Processed'

%% Import Abs data
% Check the unit of the abs data and change the wavelength in the
% fdomcorrect.m
cd ('F:\NTU\Research\Sarawak FDOM analysis_new\1cm\Abs(per meter)')
filetype = 'Abs'; ext = 'csv'; RangeIn = 'A1..B671'; display_opt = 0; outdat = 1;
[S_abs,W_abs,wave_abs,filelist_abs]=readinscans(filetype,ext,RangeIn,display_opt,outdat);

%% Import lab reagent blank data (for blank subtraction)
cd('F:\NTU\Research\Sarawak EEM\Processed\March\Lab Reagent Blank')
filetype=1;ext = 'csv';RangeIn='A1..AR157';headers=[1 1];display_opt=0;outdat=1;
[X_lrb,Emmat_lrb,Exmat_lrb,filelist_lrb,outdata_lrb]=readineems(filetype,ext,RangeIn,headers,display_opt,outdat);
Ex_lrb=Exmat_lrb(1,:); %Since all files have the same excitation wavelengths
Em_lrb=Emmat_lrb(:,1); %Since all files have the same emission wavelengths

 %% Import Di Water blank scan (for March water raman unit calibration)
% cd('F:\NTU\Research\Sarawak EEM\Processed\March\Di Blank for Raman Unit')
% filetype=1;ext = 'csv';RangeIn='A1..AR157';headers=[1 1];display_opt=0;outdat=1;
% [X_wr,Emmat_wr,Exmat_wr,filelist_wr,outdata_wr]=readineems(filetype,ext,RangeIn,headers,display_opt,outdat);
% Ex_wr=Exmat_wr(1,:); %Since all files have the same excitation wavelengths
% Em_wr=Emmat_wr(:,1); %Since all files have the same emission wavelengths

%% Import water Raman scan

cd('F:\NTU\Research\Photodegradation experiment\Maludam photodegradation May-2018\Maludam phbl FDOM June 2018\Maludam photodegradation experiment\Water Raman scan')
filetype = 'UV350';
ext = 'csv';
RangeIn = 'A2..B87';
display_opt=0;outdat=1;
[S_wr,W_wr,wave_wr,filelist_wr]=readinscans(filetype,ext,RangeIn,display_opt,outdat);  %creat a matrix of Raman Scan called S_R
RamEx=350; %note landa <> 350, so we will need to calculate RamOpt below


%% Import EEM files
cd('F:\NTU\Research\Sarawak EEM\Processed\March\FDOM EEM')   %make the EEMs subfolder in demo your current directory

filetype=1; ext = 'csv'; RangeIn='A1..AR157'; headers=[1 1];  %ext is format, like 'xls','cvs'
display_opt=0;  %to speed up, do not display each EEM
outdat=1;      %to save data
[X,Emmat,Exmat,filelist_eem,outdata]=readineems(filetype,ext,RangeIn,headers,display_opt,outdat);  %readineems
% readineems function can read in each excel file of the EEM and then
% assemble them together and generate a 3-D matrix, so here the two
% dimensions are excitation wavelength and emission wavelength and the
% third dimension is the samples: (SAMPLE*Emission*Excitation)
Ex = Exmat (1,:);
Em = Emmat (:,1);

%% Sample Log
cd('F:\NTU\Research\Sarawak EEM\Processed\March')
[LogNUM,LogTXT]=xlsread('SampleLog_March.xlsx',1);     %xlsread, LogNum can only read the numeric data, LogText can only read the text
                                                                     % 'SF-6-P is the specified spreadsheet you want in the excle file'
%Get filenames and text data from LogTXT
%Log_Date=LogTXT(:,2);
%Log_Cruise=LogTXT(:,3);
%Log_Site=LogTXT(:,4);
Log_EEMfile=LogTXT(:,5);
Log_ABSfile=LogTXT(:,1);
Log_lrbfile=LogTXT(:,2);
Log_wrfile=LogTXT(:,6);
Log_SampleID=LogTXT(:,7);
%Log_RAMfile=LogTXT(:,10);
%Log_EmCorfile=LogTXT(:,11);
%Log_ExCorfile=LogTXT(:,12);
%Log_Qw=LogTXT(:,14);

%Get numeric data from LogNUM
%RepNo=LogNUM(:,5);
%SampNo=LogNUM(:,6);
%QSdata=LogNUM(:,13);
%df=LogNUM(:,15);       %df = dilution factor
% Log_QSslope=cellstr(char('Qsslope',num2str(QSdata)));
% Log_SampNo=cellstr(char('SampNo',num2str(SampNo)));
% Log_Replicates=cellstr(char('Replicate',num2str(RepNo)));
% Log_DilFac=cellstr(char('DF',num2str(RepNo)));

%% Reconcile various datasets
% EEMs are paired with their corresponding absorbance scans and Raman scans
Pair_EEM_log=[Log_EEMfile Log_EEMfile]; %used for all numeric information in the log
Pair_EEM_wr=[Log_EEMfile Log_wrfile];
Pair_EEM_abs=[Log_EEMfile Log_ABSfile];
Pair_EEM_lrb=[Log_EEMfile Log_lrbfile];
%Pair_EEM_qsuv=[Log_EEMfile Log_Qw];

%% Obtain matching datasets - from loaded datasets
cd('F:\NTU\Research\Sarawak EEM\Processed\March\Corrected EEMs')
[Sabs Sabs_filenames]=matchsamples(filelist_eem,filelist_abs,Pair_EEM_abs,X,S_abs);% ABS scans that match filelist_eem
Swaterraman=matchsamples(filelist_eem,filelist_wr,Pair_EEM_wr,X,S_wr);        % Raman scans that match filelist_eem
Slabreagentblk=matchsamples(filelist_eem,filelist_lrb,Pair_EEM_lrb,X,X_lrb); % matching blank EEMs
SampleID_March=matchsamples(filelist_eem,Log_EEMfile(2:end,:),Pair_EEM_log,X,Log_SampleID(2:end,:)); %Match sample ID
%Sqsuv275=matchsamples(filelist_eem,filelist_qsuv275,Pair_EEM_qsuv,X,S_qsuv275); % QS blanks at Ex=275nm
%Sqsuv350=matchsamples(filelist_eem,filelist_qsuv350,Pair_EEM_qsuv,X,S_qsuv350); % QS blanks at Ex=350nm

% Obtain matching datasets - from the Sample log

% numbers only                                          
%Q=matchsamples(filelist_eem,Log_EEMfile(2:end,:),Pair_EEM_log,X,QSdata);              % QS slopes at 350/450
%sampleID=matchsamples(filelist_eem,Log_EEMfile(2:end,:),Pair_EEM_log,X,SampNo);       % SampNo matching filelist_eem
%replicates=matchsamples(filelist_eem,Log_EEMfile(2:end,:),Pair_EEM_log,X,RepNo);      % RepNo matching filelist_eem
%dilfac=matchsamples(filelist_eem,Log_EEMfile(2:end,:),Pair_EEM_log,X,df);             % df matching filelist_eem

% text only - need to remove headers, e.g. Log_Site(2:end,:);  
%sites=matchsamples(filelist_eem,Log_EEMfile(2:end,:),Pair_EEM_log,X,Log_Site(2:end,:));     % sites matching filelist_eem
%cruises=matchsamples(filelist_eem,Log_EEMfile(2:end,:),Pair_EEM_log,X,Log_Cruise(2:end,:)); % cruises matching filelist_eem
%dates=matchsamples(filelist_eem,Log_EEMfile(2:end,:),Pair_EEM_log,X,Log_Date(2:end,:));     % dates matching filelist_eem


%% FDOM correction
cd('F:\NTU\Research\Sarawak EEM\Processed\March\Corrected EEMs')

W = vertcat(wave_wr,Swaterraman)  % Extract the water raman scan at 350nm from the water blank
A = vertcat (wave_abs,Sabs);   %Preparing the matrix A for inner filter effect correction using Absorbance method
RamOpt = [350 371 428];      %use the default Raman Integration range
Emcor=horzcat(Em,ones(156,1));
Excor=horzcat(Ex',ones(43,1));  %If you want to do the fdomcorrection, there has to be a Emcor and a Excor
[XcRU2, Arp2, IFCmat2, BcRU2, XcQS2, QS_RU2]=fdomcorrect(X,Ex,Em,Emcor,Excor,W,RamOpt,A,Slabreagentblk,[],[],[]);


%% As for diluted samples, multiply 10
%Marchfile=cellstr(filelist_eem_March);
for i=1:7
    XcRU2(i,:,:)=XcRU2(i,:,:)*10;       %Maludam
end

XcRU2(23,:,:)=XcRU2(23,:,:)*10;
XcRU2(25,:,:)=XcRU2(25,:,:)*10;
XcRU2(26,:,:)=XcRU2(26,:,:)*10; 

for i=33:34
    XcRU2(i,:,:)=XcRU2(i,:,:)*10;
end

for i=38:43
    XcRU2(i,:,:)=XcRU2(i,:,:)*10;
end

for i=50:54
    XcRU2(i,:,:)=XcRU2(i,:,:)*10;
end

XcRU2(56,:,:)=XcRU2(56,:,:)*10;
%% Finally save the data
XcRU2_March=XcRU2;
filelist_eem_March=filelist_eem;

cd('F:\NTU\Research\Sarawak EEM\Processed\March\Corrected EEMs')
save CorrectedMarchDataset.mat XcRU2_March Ex Em filelist_eem_March SampleID_March
xlswrite('OutData_FDOM.xls',cellstr(filelist_eem),'filelist')
