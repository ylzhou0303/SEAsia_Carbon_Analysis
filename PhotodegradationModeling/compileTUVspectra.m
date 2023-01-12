%% in this file, I compile the TUV outputs of the hourly resolved spectra into a 3D dataset

cd('F:\NTU\Research\Photodegradation experiment\TUV model\TUV outputs\Singapore cloud corrected');
filelist = dir('F:\NTU\Research\Photodegradation experiment\TUV model\TUV outputs\Singapore cloud corrected\*.txt');

%%
AllSpectra = zeros(570,6,288)*nan;

SZA = table();
SZA.sza = ones(288,1) * nan; %save the solar zenith angle

for k = 1:288
    rawfile = fopen(filelist(k).name); %open the txt file
    outfile = fopen('joey.txt','w'); %create a file to store the data part of the txt file
    
    while ~feof(rawfile) %if it's not the end of the file
        tline = fgetl(rawfile);   %read the line from the raw file
        if ~isempty(tline) && size(tline,2)>2 %if this is not an empty line
            
            if double(tline(2)) >= 48 && double(tline(2)) <= 57  %if this is number, then save this line to joey
                fprintf(outfile, '%s\n\n',tline);
            elseif strcmp(tline(2),'s')   %if this line starts with "s" (solar zenith angle), save the value
                tempstr = tline(27:34);
                SZA.sza(k) = str2double(tempstr);
                SZA.UTC(k) = cellstr(filelist(k).name(1:10));
            end
        end
    end
    
    fclose(rawfile);
    fclose(outfile);
    AllSpectra(:,:,k) = importdata('joey.txt');
    
end


%%
AllSpectra_SG_cloud = AllSpectra;
SZA_SG = SZA;
%%
save('F:\NTU\Research\Photodegradation experiment\TUV model\TUV outputs\AllSpectra_SG_cloud.mat','AllSpectra_SG_cloud','SZA_SG');

%% plot the diurnal cycle
figure;hold on
for i = 25:48
    plot(AllSpectra(:,1,i),AllSpectra(:,6,i));
    pause(1)
end


%% plot seasonal cycle
figure;hold on
for i = 5:24:288
    plot(AllSpectra(:,1,i), AllSpectra(:,6,i));
end

xlim([200 900])
xlabel('wavelength(nm)')
ylabel('solar irradiance(W m^-^2 nm^-^1)')
title('SG cloud corrected')
legend('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')