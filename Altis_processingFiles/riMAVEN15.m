function riMAVEN15
%get the MAVEN data and the sample information squared away
%KL 10/20/2016; KL 11/14/2017 add in the ability to merge based on names
%KL 12/11/2018 - working on year 3 of BIOS-SCOPE data
%KL 2/21/2020 working on Pro8, define SRM files here 
%KL 3/16/2021 updating to use Melissa's new names
%KL 4/10/2021 - updated to use three runs on TSQ
%KL 7/5/2021 - change three standard curves (pos/neg/Mix8), but two SRM lists
%%this version will gather needed MATLAB files into the _riMAVENfiles 
%%folder so that I can have multiple functions without making one
%%giant file that has everything (considerMAVEN is getting unwieldy)
%%fxns in the folder are: catCSV, considerMAVEN_stepA/B/C, useErrors,getErrors
%KL 1/11/2023 - with the Altis we can do polarity switching and combine pos
%and neg mode, however, we still have mix 8 separate from mix1-7 and this
%requires some code shuffling.
%KL 1/4/2024 updating to remove one contaminated sample
%KL 3/12/2024 - adding sample files to show expected format.
clear
addpath('_riMAVENfiles','-end') %add the folder with the various functions to the path
%housecleaning - when testing this is leaving files behind that are annoying
fclose('all');
delete('tempFile_mix8.csv');
delete('tempFile_both.csv');

%%%set the file names up front:
NameOfFile = 'someDataFile.2024.03.12_matrix.mat';

%which curve to use? MQ or matrix
useCurve = 'matrix';

%sequence file from the TSQ; baseDir is *not* the one with sequence, one
%level higher; 
baseDir = 'sampleFiles'; 

%wDir = strcat(baseDir,filesep,'sequence_fromMethods');
wDir = baseDir; %note that KL usually had things in different folders, use this for GitHub
fName = 'mtab_sampleSequenceAltis.2024.xlsx';
sampleInfoFile = [wDir filesep fName];
clear wDir fName

%setup some variables for the TSQ processing:
setQuality = 0.1;
warning('off', 'MATLAB:codetools:ModifiedVarnames');

%ignore these labeled standards (will be in data, but not as [], peak areas)
labeledCpds = {'cholic acid d4 mix8','cholic acid d4 neg',...
    'd2-biotin neg mix8','d2-biotin pos mix8',...
    'd2-biotin neg','d2-biotin pos',...
    'glutamic acid d3 mix8 neg','glutamic acid d3 mix8 pos',...
    'glutamic acid d3 neg','glutamic acid d3 pos',...
    'indole-3-acetic acid d7 mix8','indole-3-acetic acid d7',...
    'succinic acid d6 mix8','succinic acid d6','succinic acid d6 neg',...
    'taurocholic acid d5 mix8','taurocholic acid d5','taurocholic acid d5 neg',...
    '4-hydroxybenzoic acid d4 neg','indole-3-acetic acid d7 pos'};

%%SRM files: three, but the Mix8 compounds are also duplicated in pos/neg
%SRM.sDir = strcat(baseDir,filesep,'SRMlist');
SRM.sDir = baseDir;
% SRM.pos = [SRM.sDir filesep 'SRMList_Pos_forElMAVEN.2021.05.11.csv']; %pos mode cpds and those in Mix8
% SRM.neg = [SRM.sDir filesep 'SRMList_Neg_forElMAVEN.2021.05.11.csv'];
SRM.mix8 = [SRM.sDir filesep 'SRMTable_mtab_mix8.2023.06.23.forElMAVEN.csv']; %different curve []
SRM.both = [SRM.sDir filesep 'SRMTable_mtab_mix1to8.2023.11.20.forElMAVEN.csv'];

%will need the metabolomics standards sheet for cleanName and extraction efficiency 
%you will need to pull this from the KujLab fileshare space as this is not
%on GitHub
fDir = 'Z:\_LabLogistics\MetabolomicsStandards';
standardFile = [fDir filesep 'StandardListCompilation_2024.02.16.xlsm'];
opts = detectImportOptions(standardFile,'Sheet','allMetabolites_sheet','NumHeaderLines',3);
MWinfo = readtable([standardFile],opts);
clear opts fDir standardFile ans

%use the new function to concatanate files
%No separate Mix8 analysis, but different curve, data are in pos/neg run
mix8.folder = strcat(baseDir, filesep, 'ElMAVENoutput');
mix8.fList = dir([mix8.folder filesep '*mix8*csv']);
mix8.tempFile = 'tempFile_mix8.csv';
catCSV(mix8.folder,{mix8.fList.name},mix8.tempFile);
%new function that will only pull the standard curves and filter data
mix8.stepOne =  considerMAVEN_stepA_v3(mix8.tempFile,sampleInfoFile,setQuality,1,'mix8',useCurve,SRM.mix8,[],[]);

% need the list of compound in Mix8
mix8cpds = setdiff(mix8.stepOne.mtabName,labeledCpds);

%positive and negative ion mode data
both.folder = strcat(baseDir, filesep, 'ElMAVENoutput');
both.fList = dir([both.folder filesep '*both*csv']);
both.tempFile = 'tempFile_both.csv';
catCSV(both.folder,{both.fList.name},both.tempFile);
both.stepOne =  considerMAVEN_stepA_v3(both.tempFile,sampleInfoFile,setQuality,1,'both',useCurve,SRM.both,mix8cpds,labeledCpds);

%%now get concentrations: both positive and negative ion mode data
[both.sNames both.kgd] = considerMAVEN_stepB_v2(both.stepOne,mix8.stepOne,labeledCpds);

clear setQuality SRM useCurve
clear baseDir 

%keep up with the housecleaning
oldData.mix8 = mix8; 
oldData.mix8cpds = mix8cpds; 
clear mix8*

%now, with polarity switching, I have all the data in one structure, no 
%merge needed
mtabNames = both.kgd.name;

tInfo = readtable(sampleInfoFile);
clear sampleInfoFile

%first, go through and iterate through the pooled samples to be sure there 
%are no duplicate names
s = strcmp(tInfo.SampleName,'BIOSSCOPE_2023_pooled');
ks = find(s==1);
for a = 1:length(ks);
    t = tInfo.SampleName(ks(a));
    tInfo.SampleName(ks(a)) = strcat('p',num2str(a),t);
    clear t
end
clear a ks a
 
%before I dive into the unknowns, remove anything that has goodData = 0
%also remove anything that contains Std Bracket
sc = contains(tInfo.sampleDetails,'Std Bracket');
sb = contains(tInfo.sampleDetails,'Blank');
k = find(tInfo.goodData==0 | sc==1 | sb==1);
tInfo(k,:) = [];
clear k sc sb

%January 4, 2024 - have one surface sample that is completely off from any
%prior values, delete it here so it does not propagate anywhere
%sample is BATS 10391, cast 14, Niskin
s = strcmp(tInfo.FileName,'mtab_BIOSSCOPE_temporal_2023_061223_81');
ks = find(s==1);
tInfo(ks,:) = [];
clear s ks

sInfo = table;
sInfo.cName = unique(tInfo.shortName);
%the first row of this will be empty, delete that
if isequal(sInfo.cName(1),{''});
    sInfo(1,:) = [];
end

%now make an empty matrix for the data...will be all numbers so no need for
%special format
mtabData = zeros(size(mtabNames,1),size(sInfo,1));

sInfo.runOrder(:,1) = 0;
sInfo.FileName(:,1) = {''};

for a = 1:size(sInfo,1);
    s = strcmp(sInfo.cName(a),tInfo.shortName);
    ks = find(s==1);
    
    %some variant of this to propagate sInfo with the basic information
    %add/edit rows as needed
    try
    sInfo.tempFlag(a) = tInfo.tempFlag(ks);
    catch
        keyboard
    end
%     sInfo.volFilt_L(a) = tInfo.volFilt_L(ks);
%     sInfo.cells_ml(a) = tInfo.cells_ml(ks);
%     sInfo.note(a) = tInfo.note(ks);
%     sInfo.extractVol_ul(a) = tInfo.extractVol_ul(ks);ri
        
    tName = tInfo.FileName(ks);
    sInfo.FileName(a,1) = tName;

    [c ia tIdx] =intersect(tName,both.sNames);

    try
    mtabData(:,a) = both.kgd.goodData(:,ks);
    catch
        keyboard
    end
    clear c ia tIdx tName
            
    clear s ks        
end
clear a
clear idx_* tInfo

%housecleaning
oldData.both = both; clear both

%pull run order from the filename
for a = 1: size(sInfo,1)
    gc = sInfo{a,'FileName'}{:}; %added {:} to deal with table output
    t = regexp(gc,'_');
    if ~isempty(t)
        sInfo.runOrder(a,1) = str2num(gc(t(end)+1:end));
    else
        sInfo.runOrder(a,1) = NaN;
    end
    clear gc t
end
clear a
 
%just do a little bit of housecleaning
if 1
    toDelete = {'1-deoxy-D-xylulose-5-phosphate',...
        'acetylserine','s-(1-2-dicarboxyethyl)glutathione pos',...
        'prostaglandin A2 neg','prostaglandin D2 neg','prostaglandin E2 neg',...
        'prostaglandin F2 neg','prostaglandin keto neg','biotin neg',...
        'tetrahydrobiopterin240','glutamic acid pos'};
    toDelete = cat(2,toDelete,labeledCpds);
elseif 0 %do nothing
    
end

[c ia ib] = intersect(toDelete,mtabNames);

mtabNames(ib,:)=[];
mtabData(ib,:)=[];
clear c ia ib
clear toDelete    

%check to make sure I am not missing anything. The new (2021) standard file
%now lists all possible SRMs so I can use that to find each compound. I
%have five columns of options though most compounds only have one SRM.
s = contains(MWinfo.Properties.VariableNames,'SRMname');
SRMmatch = find(s==1);
allSRMs = MWinfo{:,SRMmatch};
clear s

allSRMs = stripName(allSRMs); %adding 10/5/2022; stripName is a function KL wrote

mtabInfo=table;
%get the row information, or send up an error if I have no match
for a = 1:length(mtabNames)
    %will not have information on labeled standards so set up a skip
    sk = ismember(labeledCpds,mtabNames(a));
    if sum(sk)==0
        s = ismember(allSRMs,stripName(mtabNames(a)));
        [r,c] = find(s==1);
        if isempty(r) % sarcosine in BIOS-SCOPE... && ~strcmp(mtabNames(a),'sarcosine')
            mtabNames(a)
            keyboard
            error('Something is wrong because I have no information on this compound')
        end

        try
            mtabInfo.SRMname(a) = MWinfo{r,SRMmatch(c)};
        catch
            %probably have duplicate SRM names after stripping the names clean, 10/15/2022
            mtabNames(a)
            keyboard
        end
        mtabInfo.cleanName(a) = MWinfo{r,'cleanName'};
        mtabInfo.KEGG(a) = MWinfo{r,'KEGG'};
        mtabInfo.CHEBI(a) = MWinfo{r,'CHEBI'};
        mtabInfo.cpdCategory(a) = MWinfo{r,'category'}; %more options here
        mtabInfo.cpdCategory2(a) = MWinfo{r,'category2'}; %more options here
        mtabInfo.MW(a) = MWinfo{r,'exactNeutralMass'};
        mtabInfo.extEff(a) = MWinfo{r,'PPL_ExtEff_percent'};

        clear s r c sk
    else
        clear sk
    end
end
clear a

%keep the names of the labeledCpds - may be useful later
oldData.labeledCpds = labeledCpds; clear labeledCpds

%housecleaning - the temp files are not getting properly removed
delete('tempFile_both.csv')
delete('tempFile_mix8.csv')

rmpath('_riMAVENfiles') %remove this...causing issues when I run multiple

%get the sample inventory and pull the volume filtered from there
%The sample inventory on aston will have the volumes filtered...and the
%nominal depths
wDir = ('Z:\_LabLogistics\Sample Inventories');
fName = '2023-12-04_BIOS-SCOPE_MasterFile.xlsx';

%only need to J bc that is where the volume filtered is 
T = readtable([wDir filesep fName],'sheet','Summary','range','A9:J10000',...
    'ReadVariableNames',true);
clear wDir fName

%for now, only need the samples that are the PPL samples
s = strcmp(T.typeOfSample,'PPL: 1gm/6ml');
ks = find(s~=1);
T(ks,:) =[];
clear s ks 
%change cast and niskin to double
T.cast = str2double(T.cast);
T.niskin = str2double(T.niskin);
T.depth = str2double(T.depth);

%setup the empty column
sInfo.volFiltered_L(1:size(sInfo,1),1) = NaN;

%%now match up the samples and the volumes filtered
for a = 1:size(sInfo,1);
    one = sInfo.cName{a};
    r = regexp(one,'_');
    cruise = one(1:r-1);
    sInfo.cruise(a) = {cruise};        %need to keep the cruise information
    if contains(one(r+1:end),'C')
        rn = regexp(one,'N');
        %match by cast and niskin
        s = strcmp(cruise,T.CruiseID);
        ks = find(s==1 & T.cast==str2double(one(r+2 : rn-1)) & T.niskin == str2double(one(rn+1:end)));
        sInfo.volFiltered_L(a) = T.volumeFiltered_L_(ks);
        %need to keep cast and niskin
        sInfo.cast(a) = str2double(one(r+2 : rn-1)) ;
        sInfo.niskin(a) = str2double(one(rn+1:end));
        clear rn s ks
    elseif contains(one,'BATS386')
        %This year someone at BIOS decided to use a new naming system for
        %some samples. We have BATS386 ....which is BATS10386
        %it's always the special cases that trip you up
        s = strcmp('10386',T.CruiseID);
        ks = find(s==1 & str2double(one(r+2:end))== T.depth);
        sInfo.volFiltered_L(a) = T.volumeFiltered_L_(ks);
        sInfo.cast(a) = T.cast(ks);
        sInfo.niskin(a) = T.niskin(ks);       
        %still want the nominal depth
        sInfo.nominalDepth(a) = str2double(one(r+2:end));
        %also change the cruise variable because that gets used later
        sInfo.cruise(a) = {'BATS10386'};
        clear s ks 
    elseif contains(one(r+1:end),'D')
        %BATS cruise, match by depth
        s = strcmp(cruise(5:end),T.CruiseID);
        ks = find(s==1 & str2double(one(r+2:end))== T.depth);
        
        sInfo.volFiltered_L(a) = T.volumeFiltered_L_(ks);
        %need to manually enter cast and niskin from BATS cruise cast sheets
        %KL 12/18/2023 - that is now done and entered in Kuj inventory
        sInfo.nominalDepth(a) = str2double(one(r+2:end));
        sInfo.cast(a) = T.cast(ks);
        sInfo.niskin(a) = T.niskin(ks);
        clear s ks
    else
        %put in a bogus number for the other samples         
        sInfo.volFiltered_L(a) = 10;        
    end
    clear one r cruise
end
clear a
        
%now do the math to convert use the volume information
%follow B21/p7 math to get concentrations in the source water :
num = mtabData*0.25; %(0.25  is the volume of the extract, in milliliters)
den = (repmat(sInfo.volFiltered_L',size(mtabNames,1),1)*1000); %1000 to convert to ml
mtabData_inSeawater_pgML = (num./den)*1000; %as pg/ml
clear num den
clear mtabData T ans        
    
   

%now use the metabolite information, first to correct for extraction efficiency
mtabData_EEcorrected_pgML = nan(size(mtabData_inSeawater_pgML));
%do not correct for things with extraction efficiency less than 2%
setMinEE = 2; %set 12/19/2023

%this is mathematically inefficient, but is easier because of the error
%remember this is still rows = metabolites; columns = samples
for a = 1:size(mtabInfo,1); %rows = metabolites
    for aa = 1:size(mtabData_inSeawater_pgML,2);    %columns = samples    
        if ~isnan(mtabInfo.extEff(a)) 
            if mtabInfo.extEff(a) >= setMinEE ; 
                %convert extEff to fraction...
                mtabData_EEcorrected_pgML(a,aa) = mtabData_inSeawater_pgML(a,aa)/(mtabInfo.extEff(a)/100);
            elseif mtabInfo.extEff(a) < setMinEE  
                %%keep these, but do not correct 
                mtabData_EEcorrected_pgML(a,aa) = mtabData_inSeawater_pgML(a,aa);
            end

        elseif isnan(mtabInfo.extEff(a));
            mtabData_EEcorrected_pgML(a,aa) = mtabData_inSeawater_pgML(a,aa);
            if aa==1;
                %no need to repeat for every sample...once per mtab is fine
                fprintf(['Need extraction efficiency for ',mtabInfo.cleanName{a},'\n']);
            end
        end
    end
    clear aa
end
clear a

%these data are already in pg/ml in **seawater** now need to go from volume to moles
mtabData_moles_EEcorrected = zeros(size(mtabData_EEcorrected_pgML));

%go through one r/c at a time, inefficient, but easier to trap errors
%remember this is still rows = metabolites; columns = samples
for r = 1:size(mtabInfo,1);
    for c = 1:size(sInfo,1);
        %stop if I have no molecular weight
        if isnan(mtabInfo.MW(r))
            error(strcat('no molecular weight for ',{' '},mtabInfo.cleanName{r}))
            keyboard
        end
        %this gives moles in seawater...will be small number
        mis = mtabData_EEcorrected_pgML(r,c)*(1/1e9)*(1/mtabInfo.MW(r));
        mtabData_moles_EEcorrected(r,c) = mis;
        clear mis   
    end
end
clear r c

%change to pM
mtabData_pM_EEcorrected = mtabData_moles_EEcorrected*1e12;
clear mtabData_moles_EEcorrected mtabData_EEcorrected_pgML mtabData_inSeawater_pgML

save(NameOfFile)
end
