function matchKujTargetedMetabolites_2023
%function matchKujTargetedMetabolites_2023
%This file will pull in data from the Altis - the 2023 data
%You will need four files to make this work
%1. One run of targeted metabolites - already done with El-MAVEN and MATLAB
%2. The discrete data file from the BIOS-SCOPE Google Drive space
%3. The prior MATLAB *mat file with Kujawinski lab metabolites: 'BIOSSCOPE_metabolites_2021.mat'
%4. The lookup table ('KujLab_metaboliteLookUpTable.xlsx') which holds the
%details on what SRM names were used in what year. 
%
%This will export a MATLAB *mat file and a CSV file that can be shared
%Krista Longnecker
% 17 November 2023; updated 31 December 2023; 5 January 2024; 26 February 2024
clear all
close all
warning('off','MATLAB:table:RowsAddedExistingVars')

%2023 data are as follow:
%% These data files are on Krista's desktop in case you need to access them
ccnDir = 'C:\Users\klongnecker\Documents\Dropbox\Current projects\Kuj_BIOSSCOPE\RawData\TSQ\';
ccnName{1} = '2023_06_BIOSSCOPE temporal run3\BIOSSCOPE_RP_Altis.2024.01.04_matrix.mat'; %where ccn = cruise, cast, niskin

%Required file #1: the 2023 data, load that in the loop that follows below

%Required file #2: what is the latest version of the discrete data file, and where is it?
%% NOTE - this file sits in the BIOS-SCOPE Google Drive, get the latest version from there
bottleFile = 'BATS_BS_COMBINED_MASTER_latest.xlsx';
wDir = 'C:\Users\klongnecker\Documents\Dropbox\Current projects\Kuj_BIOSSCOPE\RawData\DataFiles_CTDandDiscreteSamples';
BottleFileInfo = [wDir filesep bottleFile];
sheetName = 'DATA'; %this sheet has the discrete data on it
%read in the Excel file
discreteData = readtable(BottleFileInfo,'sheet',sheetName);
%tidying up
clear bottleFile wDir BottleFileInfo sheetName

%Required filed #3: the MAT file with the 2016 to 2021 metabolites; metabolite concentrations are in pM
%% this file is now in the /dataFiles folder on GitHub, you may need a new path
eDir = '..\dataFiles';
eFile = 'BIOSSCOPE_metabolites_2021.mat';
existingMAT = load([eDir filesep eFile]); %
clear eDir eFile 

%Required file #4: the look up table (as Excel file) with the details on SRM names used in
%each year of the project. Names have changed over the years, new
%compounds have been added, with new instruments detection limits also
%change (note this may be a duplicate from what will end up in
%existingMAT.LUtable, but keep here in case the LU table has changed.
%% this file is now in the /dataFiles folder on GitHub, you may need a new path
tDir = '..\dataFiles';
tFile = 'KujLab_metaboliteLookUpTable.xlsx';
LUtable = readtable([tDir filesep tFile]);
clear tDir tFile 

%put in some checks. Frist, if the two versions of LUtable are different,
%there is an issue
if ~isequal(size(LUtable),size(existingMAT.LUtable))
    error('Note from Krista: you cannot proceed because seem to have new metabolites')
end
%also make sure the discrete data size matches
if ~isequal(size(discreteData),size(existingMAT.discreteData))
    error('Note from Krista: you cannot proceed because have more discrete data')
end

%housecleaning and pull out the parts I want from existingMAT
mtabData_pM_EEcorrected = existingMAT.mtabData_pM_EEcorrected;
clear existingMAT


%only one batch in 2023, keep loop in case it's useful in the future
for a2023 = 1
    %Required file #1: Go get the *mat file that resulted from one run on
    %the Altis
    %will be metabolomics data - with cruise/cast/niskin 
    oneRun = load(strcat(ccnDir,filesep,ccnName{a2023})); 
    
    %%do some housecleaning on the TSQ data before proceeding
    %get rid of the oddballs in the TSQ data
    %%do some housecleaning on the TSQ data before proceeding
    %get rid of the pooled samples and the monster samples
    sp = contains(oneRun.sInfo.cName,'pooled');
    sm = contains(oneRun.sInfo.cName,'Monster');
    sb = contains(oneRun.sInfo.cName,'MQ');
    sf = contains(oneRun.sInfo.cName,'Flow');
    gm = find(sp==1 | sm==1 | sb==1 | sf==1);
    oneRun.sInfo(gm,:)=[];
    oneRun.mtabData_pM_EEcorrected(:,gm) =[];
    clear gm sp sm sb sf
    
%     %also, for now, delete the one BATS sampling where I have no cast sheet
%     si = contains(oneRun.sInfo.cName,'BATS10388');
%     gm = find(si==1);
%     oneRun.sInfo(gm,:)=[];
%     oneRun.mtabData_pM_EEcorrected(:,gm) =[];
%     clear gm si
    
    
    %strip out the standards
    toDelete = {'4-hydroxybenzoic acid d4','cholic acid d4',...
        'd2-biotin neg','d2-biotin pos',...
        'glutamic acid d3 neg','glutamic acid d3 pos',...
        'indole-3-acetic acid d7','succinic acid d6','taurocholic acid d5'};
    [c ia ib] = intersect(stripName(toDelete),stripName(oneRun.mtabNames));

    oneRun.mtabData_pM_EEcorrected(ib,:) = [];
    oneRun.mtabNames(ib)=[];
    clear c ia ib toDelete
      
    %the sample labeling mismatch is an issue - the samples the Kujawinski lab
    %gets for the BATS cruise are labeled this way: BATS##### with depth in
    %meters. The discrete datafile has cruise as the ship information (which Kuj
    %lab does not see). Hence, the first bit is to make a dummy variable that 
    %is just the five digit cruise ID
    tempList = discreteData.ID;
    for a = 1:size(tempList,1)
        one = tempList{a};
        forKujMerge{a,1} = one(1:5);
        clear one
    end
    clear a tempList

    %go through one sample at a time, and one metabolite at a time.
    %This is inefficient, but easier to check the indexing.
    %The Kuj data is metabolites in rows and samples in columns; 
    %the discrete data for BIOS-SCOPE has samples in rows, so put the
    %metabolites in columns.

    for a = 1:size(oneRun.sInfo,1)
            oneCruise = oneRun.sInfo.cruise(a);
            oneCast = oneRun.sInfo.cast(a);
            oneNiskin = oneRun.sInfo.niskin(a);
            %now look for these three things in the discrete data file
            s = strcmp(oneCruise,discreteData.Cruise_ID);
            step = char(oneCruise);
            sm = strcmp(step(5:end),forKujMerge);
            clear step 

            ks = find(s==1 & discreteData.Cast == oneCast & discreteData.Niskin == oneNiskin);
            ksm = find(sm==1 & discreteData.Cast == oneCast & discreteData.Niskin == oneNiskin);        
        for aa = 1:size(LUtable,1);
            %now figure out which metabolite I am working on 
            oneMetabolite = LUtable.SRMname_2023(aa);
            %that metabolite name is in the lookup table, find it in oneRun
            so = strcmp(stripName(oneMetabolite),stripName(oneRun.mtabNames));
            kso = find(so==1);
                %can be missing from oneRun when it only appeared in one of
                %the batches
            if ~isempty(oneMetabolite{:}) && ~isempty(kso)
                if ~isempty(ks)
                    %this is the easiest case - cruise will match directly
                    mtabData_pM_EEcorrected(ks,aa) = oneRun.mtabData_pM_EEcorrected(kso,a);
                elseif ~isempty(ksm)
                    %for the BATS cruises, use forKujMerge
                    mtabData_pM_EEcorrected(ksm,aa) = oneRun.mtabData_pM_EEcorrected(kso,a); 
                elseif isequal(oneCruise,{'BATS10380'})
                    %vials labeled 10380, collected on 20380 (per conversation 
                    %with KL and Rachel Parsons); and both details will
                    %be missing from the BIOS-SCOPE discrete file. I think
                    %there are too many unknowns, let's skip these                    
                else
                    warning('something is wrong, did not get a match for cruise/cast/niskin')
                    keyboard
                end
            else
                %do nothing, do not have this metabolite this year
            end
            clear oneMetabolite so kso
        end
        clear aa
        clear s sm ks ksm oneCruise oneCast oneNiskin
    end
    clear a

    %tidying up (set the following to 1 if all the data check out)
    if 1
        clear forKujMerge
    end
      

end
clear a2023
clear ccnDir ccnName

%save the MATLAB file, the person who follows after me will need some MATLAB skills.
save('..\dataFiles\BIOSSCOPE_metabolites_2023.mat');

%this is the most recent file, so I will also export an Excel file that the others
%in the BIOS-SCOPE project can read. I am using Excel so I can export 
%the metabolite information AND the metabolite data in one file
%now add in the clean metabolite names so people know what metabolite is what
eFilename = '..\KujawinskiWHOI_targetedMetabolites.2024.02.26.xlsx';

forExport_New_ID = array2table(discreteData.New_ID); %update to use New_ID February 2024

%BIOS-=SCOPE uses -999 for missing data, so change the NaNs to -999
i = isnan(mtabData_pM_EEcorrected);
mtabData_pM_EEcorrected(i) = -999;
clear i
forExport_data = array2table(mtabData_pM_EEcorrected);
forExport_data.Properties.VariableNames = LUtable.cleanName;

forExport = cat(2,forExport_New_ID,forExport_data);
if isequal(forExport.Properties.VariableNames{1},'Var1') %lazy hack to change name
    forExport.Properties.VariableNames(1) = {'New_ID'};
end

writetable(forExport,eFilename,'Sheet','metaboliteData') %this is the data
writetable(LUtable,eFilename,'Sheet','metaboliteInfo'); %also export the mtabInfo
    
end

%need this function, just make it an internal function
function newOne = stripName(one)
%function newOne = stripName(one)
%the metabolite names are really tedious, strip all the annoying pieces 
%so that the matching is easier. 
%Krista Longnecker 10/5/2022; the BC lists have yet another special character 3/1/2023
%note - the easist way to find new oddball characters is to get a string of
%char with the special character and do this:
%double(str) -- there will be one number per char, and you can use that
%below in the first if statement to delete the special character
%There are three steps needed to make these names manageable
%1. strip out all the apostrophes --> char(39)
%1b. strip out the hyphens --> char(45)
%1c. strip out the apostrophe --> char(8242)
%2. make names all lower
%3. deblank the names
    if ischar(one)  
        %char will only be one metabolite
        %do nothing...already ready for next step
        one(one==char(39))=[];
        one(one==char(45))=[];
        %adding more...I wonder if there is a way to remove all special characters.
        one(one==char(8242))=[]; %3/1/2023 - this is in the BC list 
        one(one==char(226))=[]; %3/1/2023 - this is in the BC list 
        one(one==char(8364))=[]; %3/1/2023 - this is in the BC list
        one(one==char(178))=[]; %3/1/2023 - this is in the BC list 
        newOne = lower(deblank((one)));    
    else
        %seems to be a cell - could be one metabolite, or a cell array of them
        %do the stripping and then give back a cell
        if length(one)==1
            one = one{:};
            one(one==char(39))=[];
            one(one==char(45))=[];
            newOne = {lower(deblank((one)))};
        else
            %have multiple names in a cell, each row is a char
            newOne = cell(size(one));
            for a = 1:size(one,1)
                for j = 1:size(one,2) 
                    t = one{a,j};
                    t(t==char(39))=[];
                    t(t==char(45))=[];
                    newOne{a,j} = lower(deblank(t));
                    clear t
                end
            end
            clear a
        end
    end

end %close the stripName function

