function matchKujTargetedMetabolites_2016to2019
%function matchKujTargetedMetabolites_2016to2019
%Match the rows for the WHOI/Kujawinski targeted metabolomics data. Start 
%code from scratch because the existing code has gotten a bit bulky over the years.
%You will need three files to make this work
%1. The combined TSQ data from 2106 until 2019 - in one *mat file
%2. The discrete data file from the BIOS-SCOPE Google Drive space
%3. The list of Kujawinski lab metabolomics standards, with extraction
% efficiency information and molecular weights
% This will export a MATLAB mat file with the results
%
% Krista Longnecker, Woods Hole Oceanographic Institution
% 17 November 2023; updated 30 December 2023; 5 January 2024
clear all
close all
warning('off','MATLAB:table:RowsAddedExistingVars')

%Required file #1: Use the existing *mat file from the Kujawinski lab to match up the
%metabolomics data - use the one file that is the merged TSQ data with cruise/cast/niskin 
%% NOTE - this file is too big for GitHub, so you will need to get it from
%% Krista's computer if you need the original file (5 January 2024)
ccnName = 'BIOS_temporal_TSQ_useSW.2022.10.15_v1.mat'; %where ccn = cruise, cast, niskin
ccnDir = 'C:\Users\klongnecker\Documents\Dropbox\Current projects\Kuj_BIOSSCOPE\RawData';
targetedMetabolites = load(strcat(ccnDir,filesep,ccnName)); 
clear ccnName ccnDir

%Required file #2: what is the latest version of the discrete data file, and where is it?
%% NOTE - this file sits in the BIOS-SCOPE Google Drive, get the latest version from there
bottleFile = 'BATS_BS_COMBINED_MASTER_2024.01.04.xlsx';
wDir = 'C:\Users\klongnecker\Documents\Dropbox\Current projects\Kuj_BIOSSCOPE\RawData\DataFiles_CTDandDiscreteSamples';
BottleFileInfo = [wDir filesep bottleFile];
sheetName = 'BATS_BS bottle file'; %this sheet has the discrete data on it
%read in the Excel file
discreteData = readtable(BottleFileInfo,'sheet',sheetName);
%tidying up
clear bottleFile wDir BottleFileInfo sheetName

%Required file #3: need the Excel file with details on the Kuj lab metabolites
%% this file is now in the /dataFiles folder on GitHub, you may need a new path
mFile = '../dataFiles';
mDir = 'StandardListCompilation_2024.01.02.xlsm';
standardFile = [mFile filesep mDir];
opts = detectImportOptions(standardFile,'Sheet','allMetabolites_sheet','NumHeaderLines',3);
MWinfo = readtable([standardFile],opts);
clear opts mFile mDir standardFile   
       
%Required file #4: the look up table (as Excel file) with the details on SRM names used in
%each year of the project. Names have changed over the years, new
%compounds have been added, with new instruments detection limits also
%change
%% this file is now in the /dataFiles folder on GitHub, you may need a new path
tDir = '../dataFiles';
tFile = 'KujLab_metaboliteLookUpTable.xlsx';
LUtable = readtable([tDir filesep tFile]);
clear tDir tFile 
    

%%do some housecleaning on the TSQ data before proceeding
%get rid of the oddballs in the TSQ data
so = strcmp(targetedMetabolites.sInfo.addedInfo,'Low');
sr = strcmp(targetedMetabolites.sInfo.addedInfo,'R2');
sm = strcmp(targetedMetabolites.sInfo.addedInfo,'MC13C20');
sb = strcmp(targetedMetabolites.sInfo.addedInfo,'BlankMQ');
ks = find(so==1 | sr == 1 | sm==1 | sb==1);
clear so sr sm sb
targetedMetabolites.sInfo(ks,:)=[];
targetedMetabolites.mtabData_inSeawater_pgML(:,ks) =[];
targetedMetabolites.mtabData_exceedLOD(:,ks) =[];

clear s ks

%also get rid of the pooled samples (and the monster sample); cruise will
%be empty for the pooled samples - so you that to filter
gm = cellfun(@isempty,targetedMetabolites.sInfo.cruise);
targetedMetabolites.sInfo(gm,:)=[];
targetedMetabolites.mtabData_inSeawater_pgML(:,gm) =[];
targetedMetabolites.mtabData_exceedLOD(:,gm) =[];
clear gm


%make choices about positive versus negative ion mode data
toDelete = {'2deoxyinosine neg','biotin neg','cytidine neg',...
    'd2phosphoglyceric acid','nacetyl muramic acid loss h2o',...
    'phosphohomoserine neg','sarcosine'};

[c ia ib] = intersect(stripName(toDelete),targetedMetabolites.mtabNames);

if ~isequal(length(c),length(toDelete))
    error('Something is wrong, these should be the same length')
end

targetedMetabolites.mtabData_inSeawater_pgML(ib,:) = [];
targetedMetabolites.mtabData_exceedLOD(ib,:) = [];
targetedMetabolites.mtabNames(ib)=[];
clear c ia ib toDelete

%Next correct for extraction efficiency, use information from the
%Excel file with all the Kuj metabolomics compounds

%I have five columns of options though most compounds only have one SRM.
s = contains(MWinfo.Properties.VariableNames,'SRMname');
SRMmatch = find(s==1);
allSRMs = MWinfo{:,SRMmatch};
clear s
allSRMs = stripName(allSRMs); %strip extra characters from the names

%match to the information in targetedMetabolites
mtabInfo = table();

%get the row information, or send up an error if I have no match
for a = 1:length(targetedMetabolites.mtabNames)
    s = ismember(allSRMs,targetedMetabolites.mtabNames(a));
    [r,c] = find(s==1);
    if isempty(r) 
        targetedMetabolites.mtabNames(a)
        error('Something is wrong because I have no information on this compound')
    end
       
    mtabInfo.SRMname(a) = MWinfo{r,SRMmatch(c)}; %if there is an error here - probably have duplicate names after stripping the name

    mtabInfo.cleanName(a) = MWinfo{r,'cleanName'};
    mtabInfo.KEGG(a) = MWinfo{r,'KEGG'};
    mtabInfo.cpdCategory(a) = MWinfo{r,'category'}; %more options here
    mtabInfo.cpdCategory2(a) = MWinfo{r,'category2'}; %more options here
    %mtabInfo.ionMode(a) =  %do not have that information there
    mtabInfo.MW(a) = MWinfo{r,'exactNeutralMass'};
    mtabInfo.extEff(a) = MWinfo{r,'PPL_ExtEff_percent'};
    
    clear s r c
end
clear a
clear allSRMs SRMmatch MWinfo

%now use the metabolite information, first to correct for extraction efficiency
targetedMetabolites.mtabData_EEcorrected_pgML = nan(size(targetedMetabolites.mtabData_inSeawater_pgML));
%do not correct for things with extraction efficiency less than 2%
setMinEE = 2; %set 12/19/2023

%this is mathematically inefficient, but is easier because of the indexing
%remember this is still rows = metabolites; columns = samples
for a = 1:size(mtabInfo,1); %rows = metabolites
    for aa = 1:size(targetedMetabolites.mtabData_inSeawater_pgML,2);    %columns = samples    
        if ~isnan(mtabInfo.extEff(a)) 
            if mtabInfo.extEff(a) >= setMinEE ; 
                %convert extEff to fraction...
                targetedMetabolites.mtabData_EEcorrected_pgML(a,aa) = targetedMetabolites.mtabData_inSeawater_pgML(a,aa)/(mtabInfo.extEff(a)/100);
            elseif mtabInfo.extEff(a) < setMinEE  
                %%keep these, but do not correct 
                targetedMetabolites.mtabData_EEcorrected_pgML(a,aa) = targetedMetabolites.mtabData_inSeawater_pgML(a,aa);
            end
            
        elseif isnan(mtabInfo.extEff(a));
            targetedMetabolites.mtabData_EEcorrected_pgML(a,aa) = targetedMetabolites.mtabData_inSeawater_pgML(a,aa);
            if aa==1;
                %no need to repeat for every sample...once per mtab is fine
                fprintf(['Need extraction efficiency for ',mtabInfo.cleanName{a},'\n']);
            end
        end
    end
    clear aa
end
clear a

%now need to go from volume to moles
targetedMetabolites.mtabData_moles_EEcorrected = zeros(size(targetedMetabolites.mtabData_EEcorrected_pgML));

%go through one r/c at a time, inefficient, but easier to trap errors
%remember this is still rows = metabolites; columns = samples
for r = 1:size(mtabInfo,1);
    for c = 1:size(targetedMetabolites.sInfo,1);
        %this gives moles in seawater...will be small number
        mis = targetedMetabolites.mtabData_EEcorrected_pgML(r,c)*(1/1e9)*(1/mtabInfo.MW(r));
        targetedMetabolites.mtabData_moles_EEcorrected(r,c) = mis;
        clear mis   
    end
end
clear r c

%change to pM
targetedMetabolites.mtabData_pM_EEcorrected = targetedMetabolites.mtabData_moles_EEcorrected*1e12;


%%now set up the match to the discrete data
%%discrete data file is rows = samples; columns = metabolites

%go through one row at a time in targetedMetabolites data. There should be
%no cases where something does not exist in the discreteData matrix.

%first, make matrix that matches the size of discreteData
mtabData_pM_EEcorrected = nan(size(discreteData,1),size(targetedMetabolites.mtabNames,1));

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

%go through one sample at a time
%The Kuj data is metabolites in rows and samples in columns; 
%the discrete data for BIOS-SCOPE has samples in rows, so put the
%metabolites in columns.
%since this is the first batch, just use all the metabolites. In future
%years will need to track which metabolite goes where
for a = 1:size(targetedMetabolites.sInfo,1)
    oneCruise = targetedMetabolites.sInfo.cruise(a);
    oneCast = targetedMetabolites.sInfo.cast(a);
    oneNiskin = targetedMetabolites.sInfo.niskin(a);
    %now look for these three things in the discrete data file
    s = strcmp(oneCruise,discreteData.Cruise_ID);
    step = char(oneCruise);
    sm = strcmp(step(5:end),forKujMerge);
    clear step 
    
    ks = find(s==1 & discreteData.Cast == oneCast & discreteData.Niskin == oneNiskin);
    ksm = find(sm==1 & discreteData.Cast == oneCast & discreteData.Niskin == oneNiskin);
    
    if ~isempty(ks)
        %this is the easiest case - cruise will match directly
        mtabData_pM_EEcorrected(ks,:) = targetedMetabolites.mtabData_pM_EEcorrected(:,a);
    elseif ~isempty(ksm)
        %for the BATS cruises, use forKujMerge
        mtabData_pM_EEcorrected(ksm,:) = targetedMetabolites.mtabData_pM_EEcorrected(:,a);
    elseif strcmp(oneCruise,'HS1319')
        %Kuj lab has samples from one Hydrostation S cruise
        s = strcmp('AE1719',discreteData.Cruise_ID);
        ks = find(s==1 & discreteData.Cast == oneCast & discreteData.Niskin == oneNiskin);
        mtabData_pM_EEcorrected(ks,:) = targetedMetabolites.mtabData_pM_EEcorrected(:,a);        
    else
        warning('something is wrong, did not get a match for cruise/cast/niskin')
        keyboard
    end
    clear s sm ks ksm oneCruise oneCast oneNiskin
end
clear a

%tidying up (set the following to 1 if all the data check out)
if 1
    clear forKujMerge
end
   
%now, expand the mtabData to match the full set of mtabs in LUtable
[c ia ib] = intersect(stripName(LUtable.SRMname_2016to2019),stripName(mtabInfo.SRMname));
existing_mtabData_pM_EEcorrected = nan(size(mtabData_pM_EEcorrected,1),size(LUtable,1));
%now populate the array with the data I have from existing
existing_mtabData_pM_EEcorrected(:,ia) = mtabData_pM_EEcorrected(:,ib);
clear c ia ib
%tidying up and rename
clear mtabData_pM_EEcorrected mtabInfo targetedMetabolites
mtabData_pM_EEcorrected = existing_mtabData_pM_EEcorrected;
clear existing_mtabData_pM_EEcorrected
    
%save the MATLAB file, the person who follows after me will need some MATLAB skills.
save('../dataFiles/BIOSSCOPE_metabolites_2016to2019.mat');

end

%need this function, make an internal function
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

