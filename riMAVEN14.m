function riMAVEN14
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
%KL 1/28/2022 stick this on GitHub bc I keep needing it in different places
clear
addpath('_riMAVENfiles','-end') %add the folder with the various functions to the path
%housecleaning - when testing this is leaving files behind that are
%annoying
fclose('all');
delete('tempFile_pos.csv');
delete('tempFile_neg.csv');
delete('tempFile_both.csv');
delete('tempFile_hil.csv');

%%%set the file names up front:
NameOfFile = 'Pro9_RP_TSQdata.2021.12.17_matrix.mat';

%whichCurve = 'matrix' or 'MQ';
useCurve = 'MQ';

%sequence file from the TSQ
baseDir = 'C:\Users\klongnecker\Documents\Current projects\Kujawinski_Chisholm_Prochlorococcus\TSQ\2021_0929_Pro9_RP';
wDir = strcat(baseDir,filesep,'sequence_fromMethods');
fName = 'mtab_Pro9_RP_2021_0929.sequenceKLcorrected_v3.xlsx';
sampleInfoFile = [wDir filesep fName];
clear wDir fName

%setup some variables for the TSQ processing:
setQuality = 0.1;
warning('off', 'MATLAB:codetools:ModifiedVarnames');

%ignore these labeled standards
labeledCpds = {'cholic acid d4 mix8','cholic acid d4',...
    'd2-biotin neg mix8','d2-biotin pos mix8',...
    'd2-biotin neg','d2-biotin pos',...
    'glutamic acid d3 mix8 neg','glutamic acid d3 mix8 pos',...
    'glutamic acid d3 neg','glutamic acid d3 pos',...
    'indole-3-acetic acid d7 mix8','indole-3-acetic acid d7',...
    'succinic acid d6 mix8','succinic acid d6',...
    'taurocholic acid d5 mix8','taurocholic acid d5',...
    '4-hydroxybenzoic acid d4'};

%%SRM files: three, but the Mix8 compounds are also duplicated in pos/neg
SRM.sDir = strcat(baseDir,filesep,'SRMlist');
SRM.pos = [SRM.sDir filesep 'SRMList_Pos_forElMAVEN.2021.05.11.csv']; %pos mode cpds and those in Mix8
SRM.neg = [SRM.sDir filesep 'SRMList_Neg_forElMAVEN.2021.05.11.csv'];
SRM.mix8 = [SRM.sDir filesep 'SRMList_Mix8_030921_UptakeIncubations_forElMAVEN_v2.csv'];

%use the new function to concatanate files
%No separate Mix8 analysis, but different curve, data are in pos/neg run
mix8.folder = strcat(baseDir);
mix8.fList = dir([mix8.folder filesep '*both*csv']);
mix8.tempFile = 'tempFile_both.csv';
catCSV(mix8.folder,{mix8.fList.name},mix8.tempFile);
%new function that will only pull the standard curves and filter data
mix8.stepOne =  considerMAVEN_stepA_v2(mix8.tempFile,sampleInfoFile,setQuality,1,'both',useCurve,SRM.mix8,[],[]);

% need the list of compound in Mix8
mix8cpds = setdiff(mix8.stepOne.mtabName,labeledCpds);

%negative ion mode data
neg.folder = strcat(baseDir);
neg.fList = dir([neg.folder filesep '*neg*csv']);
neg.tempFile = 'tempFile_neg.csv';
catCSV(neg.folder,{neg.fList.name},neg.tempFile);
neg.stepOne =  considerMAVEN_stepA_v2(neg.tempFile,sampleInfoFile,setQuality,1,'negative',useCurve,SRM.neg,mix8cpds,labeledCpds);

%positive ion mode data
pos.folder = strcat(baseDir);
pos.fList = dir([pos.folder filesep '*pos*csv']);
pos.tempFile = 'tempFile_pos.csv';
catCSV(pos.folder,{pos.fList.name},pos.tempFile);
pos.stepOne =  considerMAVEN_stepA_v2(pos.tempFile,sampleInfoFile,setQuality,1,'positive',useCurve,SRM.pos,mix8cpds,labeledCpds);
clear useCurve

%%now get concentrations: positive ion mode data
[pos.sNames pos.kgd] = considerMAVEN_stepB_v2(pos.stepOne,mix8.stepOne,labeledCpds);

%%now get concentrations: negative ion mode data
[neg.sNames neg.kgd] = considerMAVEN_stepB_v2(neg.stepOne,mix8.stepOne,labeledCpds);

clear setQuality SRM

%keep up with the housecleaning
oldData.mix8 = mix8; 
oldData.mix8cpds = mix8cpds; 
clear mix8*

%now, I have two structures with data (pos, neg) - merge these together
%(can use the old code for this)

%get all the metabolite names...should be a unique list even when I merge
%across positive and negative ion mode, but check that
mtabNames = sort(cat(1,neg.kgd.name,pos.kgd.name));
if length(unique(mtabNames)) ~= length(mtabNames)
    error('Something is wrong - duplicate names in the list of metabolites')
end

tInfo = readtable(sampleInfoFile);
clear sampleInfoFile

%first, go through and iterate through the pooled samples to be sure there 
%are no duplicate names
s = strcmp(tInfo.SampleName,'Pro9_RP_pool pos');
ks = find(s==1);
for a = 1:length(ks);
    t = tInfo.SampleName(ks(a));
    tInfo.SampleName(ks(a)) = strcat('p',num2str(a),t);
    clear t
end
clear a ks a
    
s = strcmp(tInfo.SampleName,'Pro9_RP_pool neg');
ks = find(s==1);
for a = 1:length(ks);
    t = tInfo.SampleName(ks(a));
    tInfo.SampleName(ks(a)) = strcat('p',num2str(a),t);
    clear t
end
clear a ks a   

%before I dive into the unknowns, remove anything that has goodData = 0
%at this point, also only want Unknowns
s = strcmp(tInfo.SampleType,'Unknown');
k = find(tInfo.goodData==0 | s ~= 1);
tInfo(k,:) = [];
clear k s

%now find the Unknowns...should have the same number for positive and
%negative ion mode bc have pruned out the different QC samples already
s = strcmp(tInfo.SampleType,'Unknown');
sp = strcmp(tInfo.ionMode,'positive');
ksp = (find(s==1 & sp==1));
sn = strcmp(tInfo.ionMode,'negative');
ksn = (find(s==1 & sn==1));

if ~isequal(length(ksp),length(ksn))
    error('Something wrong, these should be the same length')
end
clear a sp sn ksp ksn


sInfo = table;
sInfo.cName = unique(tInfo.shortName);
%the first row of this will be empty, delete that
if isequal(sInfo.cName(1),{''});
    sInfo(1,:) = [];
end

%now make an empty matrix for the data...will be all numbers so no need for
%special format
mtabData = zeros(size(mtabNames,1),size(sInfo,1));
%need to track some additional details:
mtabDetails = table();

%get the index for rows for positive AND negative mtabs:
[c idx_posNew idx_posOld] = intersect(mtabNames,pos.kgd.name);
[c idx_negNew idx_negOld] = intersect(mtabNames,neg.kgd.name);

mtabDetails.mode(idx_posNew,1) = {'positive'};
mtabDetails.mode(idx_negNew,1) = {'negative'};

sInfo.runOrder_pos(:,1) = 0;
sInfo.runOrder_neg(:,1) = 0;

sInfo.FileName_pos(:,1) = {''};
sInfo.FileName_neg(:,1) = {''};

for a = 1:size(sInfo,1);
    s = strcmp(sInfo.cName(a),tInfo.shortName);
    ks = find(s==1);
    if length(ks) ~= 2
        sInfo.cName(a)
        error('Something is wrong, should be two of each')
    end
    
    %some variant of this:
    for aa = 1:2
        %propagate sInfo with the basic information, only do once
        %no need to completely fill this in bc not working on data analysis
        if 1 %minimal added information for Cara
            if aa == 1            
                sInfo.sample(a) = tInfo.sample(ks(aa));
                sInfo.strain(a) = tInfo.strain(ks(aa));
                sInfo.media(a) = tInfo.media(ks(aa));
                sInfo.replicate(a) = tInfo.replicate(ks(aa));
                % sInfo.fluid(a) = tInfo.fluid(ks(aa));
                % sInfo.volSent_ml(a) = tInfo.volSent_ml(ks(aa));
                % sInfo.cellsML(a) = tInfo.cellsML(ks(aa));
                % sInfo.cellsFiltered(a) = tInfo.cellsFiltered(ks(aa));
            end
        end
        
        im = tInfo.ionMode{ks(aa)};
        if isequal(im,'positive')
            tName = tInfo.FileName(ks(aa));
            sInfo.FileName_pos(a,1) = tName;

            [c ia tIdx] =intersect(tName,pos.sNames);
            mtabData(idx_posNew,a) = pos.kgd.goodData(idx_posOld,tIdx);
            clear c ia tIdx tName
            
        elseif isequal(im,'negative')
            tName = tInfo.FileName(ks(aa));
            sInfo.FileName_neg(a,1) = tName;

            [c ia tIdx] =intersect(tName,neg.sNames);
            mtabData(idx_negNew,a) = neg.kgd.goodData(idx_negOld,tIdx);
            clear c ia tIdx tName
        else 
            error('Something wrong')
        end
        clear im
    end
    clear aa s ks        
end
clear a

clear idx_* tInfo

%housecleaning
oldData.pos = pos; clear pos
oldData.neg = neg; clear neg

for a = 1: size(sInfo,1)
    %do positive ion mode first
    gc = sInfo{a,'FileName_pos'}{:}; %added {:} to deal with table output
    t = regexp(gc,'_');
    if ~isempty(t)
        sInfo.runOrder_pos(a,1) = str2num(gc(t(end)+1:end));
    else
        sInfo.runOrder_pos(a,1) = NaN;
    end
    clear gc t
    
    %then negative ion mode first
    gc = sInfo{a,'FileName_neg'}{:}; %added {:} to deal with table output
    t = regexp(gc,'_');
    if ~isempty(t)
        sInfo.runOrder_neg(a,1) = str2num(gc(t(end)+1:end));
    else
        sInfo.runOrder_neg(a,1) = NaN;
    end
    clear gc t
end
clear a
 
if 0
    % %making decisions about pos/neg mtabs, and pruning compounds I do not
    % trust (sucrose/trehalose), or want (prostaglandins)
    toDelete = {'NAD neg','adenine neg','adenosine 5''-monophosphate neg'...
        'desthiobiotin neg','guanosine neg',...
        'hemin I','n-acetyl muramic acid neg','pantothenic acid pos',...
        's-(5''-adenosyl)-L-homocysteine neg', 'xanthine neg',...
        'xanthosine pos','taurocholic acid d5',...
        'glutamic acid d3 pos','cytidine neg',...
        'sn-glycerol 3-phosphate pos ','chitobiose pos',...
        'biotin neg','folic acid neg','meso26diaminopimelicacid',...
        '2deoxyinosine neg','aspartic acid pos',...
        'sucrose_341 neg','sucrose_387 neg','trehalose_341 neg',...
        'trehalose_387 neg','Sucrose pos','Trehalose pos',...
        'prostaglandin A2 neg','prostaglandin D2 neg','prostaglandin E2 neg',...
        'prostaglandin F2 neg','prostaglandin keto neg','D-2-phosphoglyceric acid',...
        'acetylserine','1-deoxy-D-xylulose-5-phosphate','meso26diaminopimelicacid',...
        '5-hydroxyLtryptophan','folynic acid pos','n-acetyl muramic acid loss h2o pos',...
        'phosphoserine pos','sn-glycerol 3-phosphate pos',...
        '2deoxyguanosine pos','3deoxyguanosine pos','fructose 6-phosphate',...
        'hydroxypyruvate','hypoxanthine','6-phosphogluconic acid','chitotriose',...
        'inosine 5''-monophosphate pos','myo-inositol neg','s-(1-2-dicarboxyethyl)glutathione pos',...
        'D-erythrose-4-phosphate','D-fructose 1_6-bisphosphate',...
        '4-hydroxybenzoic acid d4','glutamic acid d3 neg','taurocholic acid d5'};  %deleted labeled standards not in Cara's data (ones for filters)};
    
    [c ia ib] = intersect(toDelete,mtabNames);
    
    mtabNames(ib,:)=[];
    mtabDetails(ib,:)=[];
    mtabData(ib,:)=[];
    clear c ia ib
    clear toDelete    
end

% 
% %go back and get the cleanName and extraction efficiency information
% %but, allow the labeled standards to be as is (in labeledCpds)
% try
%     fDir = 'Z:\_LabLogistics\MetabolomicsStandards';
%     standardFile = [fDir filesep 'StandardListCompilation_2021.09.06.xlsm'];
%     ImportOptions_ = detectImportOptions([standardFile], 'NumHeaderLines', 4,...
%         'Sheet','allMetabolites_sheet','VariableNamesRange','5:5');
% catch %if I don't have VPN on, set up to use file on my desktop
%     fDir = 'C:\Users\klongnecker\Documents\Current projects\MSdataAnalysis\TSQ_standardsList';
%     standardFile = [fDir filesep 'StandardListCompilation_2021.09.06.xlsm'];
%     ImportOptions_ = detectImportOptions([standardFile], 'NumHeaderLines', 4,...
%         'Sheet','allMetabolites_sheet','VariableNamesRange','5:5');
%     fprintf('Careful - this is the metabolomics list from a desktop: may not be current\n')
% end
% MWinfo = readtable([standardFile],ImportOptions_); clear ImportOptions_
% clear opts
% 
% s = contains(MWinfo.Properties.VariableNames,'SRMname');
% SRMmatch = find(s==1);
% allSRMs = MWinfo{:,SRMmatch};
% clear s 
% 
% for a = 1:size(mtabNames,1)
%     if contains(mtabNames{a},labeledCpds)
%         %leave as is
%         mtabDetails.SRMname(a) = mtabNames(a);
%         mtabDetails.cleanName(a) = mtabNames(a);
%         mtabDetails.extractionEfficiency(a) = NaN;
%     else
%         s = ismember(lower(allSRMs),lower(mtabNames(a)));
%         [r,c] = find(s==1);
%         if isempty(r) % sarcosine in BIOS-SCOPE... && ~strcmp(mtabNames(a),'sarcosine')
%             mtabNames(a)
%             error('Something is wrong because I have no information on these compounds')
%         end
%         mtabDetails.SRMname(a) = MWinfo{r,SRMmatch(c)};
%         mtabDetails.cleanName(a) = MWinfo{r,'cleanName'};
%         mtabDetails.extractionEfficiency(a) = MWinfo{r,'PPL_ExtEff_percent'};
%         clear r c s
%     end
% end
% clear a MWinfo fDir standardFile allSRMs SRMmatch

oldData.labeledCpds = labeledCpds; clear labeledCpds

%for Pro filters, do this:
%for glutathione oxidized, need to set the filter samples to zeros
fprintf('Glutathione oxidized in filters set to zero\n')
sr = strcmp(mtabNames,'glutathione oxidized pos'); %update name 12/2021
ksr = find(sr==1);
sc = strcmp(sInfo.sample,'filter');
ksc = find(sc==1);

mtabData(ksr,ksc) = 0;
clear sr ksr sc ksc
% 
% pos_notCurrent = pos;
% clear pos

%housecleaning - the temp files are not getting properly removed
delete('tempFile_both.csv')
delete('tempFile_pos.csv')
delete('tempFile_neg.csv')
save(NameOfFile)
end