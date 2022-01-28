function [exportCurves] = considerMAVEN_stepA_v2(CSVfileFromMAVEN,sampleInfoFile,setQuality,pruneData,ionMode,curveType,SRMfile,mix8cpds,labeledStandards)
%function [exportCurves] = considerMAVEN_stepA_v2(CSVfileFromMAVEN,sampleInfoFile,setQuality,pruneData,ionMode,curveType,SRMfile,mix8cpds,labeledStandards)
%this is based on considerMAVEN, but will (1) get the standard curves for 
%the compounds, and (2) uses the steps to filter good/bad data from MAVEN
%curves which will then be applied to the samples which were analyzed in
%the positive and negative ion mode run.
%input files:
%1: CSVfileFromMAVEN: the name of the CSV file from MAVEN
%2: sampleInfoFile: the details on where the file names for the samples, 
%blanks, and standards 
%3: setQuality: the required level for a peak in the standard curve to be
%considered as good. Have been using 0.1, but might want to be more lenient
%at times
%working on using the detailed information from MAVEN to construct
%calibration curves. Want to set specific criteria before I will allow a
%peak to be used in a calibration curve. These criteria will be: (1) the
%quality score from MAVEN, (2) the signal:noise ratio of the peak in MAVEN,
%and (3) something that considers the range of values found in the samples
%4. pruneData, set to 1 if I want to only consider files marked as goodData
%5. ionMode is positive/negative/both --> three runs on the TSQ
%6. SRMfile is the name of the CSV file with the SRM list, one per ion mode
%in the sample info file
%7. mix8cpds are the compounds in mix8, don't calculate curves for those IF
%the data are pos or neg mode
% This code was previously part of considerMAVEN
%KL 6/30/2021

% warning('off','MATLAB:table:ModifiedAndSavedVarnames')
% warning('off','stats:dataset:ModifiedVarnames')
warning('off')

switch ionMode
    case 'positive'
        %maybe this will be easier to read:
        special = table();
        special.cases(1) = {'2-deoxyinosine_Na pos'}; special.match(1) = {''};
        special.cases(2) = {'asparagine'}; special.match(2) = {'asparagine confirm 2'};
        special.cases(3) = {'isoleucine'}; special.match(3) = {'isoleucine confirm 2'};
        
    case 'negative'
        special = table();
        special.cases(1) = {'desthiobiotin neg'}; special.match(1) = {''};
        special.cases(2) = {'3-methyl-2-oxobutanoic acid'}; special.match(2) = {''};
        special.cases(3) = {'4-methyl-2-oxopentanoic acid'}; special.match(3) = {''};
        special.cases(4) = {'3-3-dimethyl-2oxobutanoic acid'}; special.match(4) = {''};
        special.cases(5) = {'2-methyl-4-oxopentanoic acid'}; special.match(5) = {''};
        special.cases(6) = {'isethionic acid '}; special.match(6) = {'isethionic acid confirm 1'}; 
        special.cases(7) = {'2oxohexanoic acid'}; special.match(7) = {''};
        
    case 'both'
        special = table();
        special.cases(1) = NaN; special.match(1) = NaN;
        
end

srm = dataset('File',SRMfile,'delimiter',',');

%easiest to find the confirm compounds
r = regexp(srm.compound,'confirm');
gm = cellfun(@isempty,r); %will be 1 when the compound is NOT labeled as confirm
k = find(gm==1);

compoundList = mat2dataset(srm.compound(k,1),'VarNames',{'name'}); 
compoundList.indexMain = 0;
compoundList.indexConfirm = 0;
clear r gm k

%go through each compound, find the index for the main peak and the confirm 
for a = 1:length(compoundList);
    %change syntax 3/16/2021
    %sn = strcmp(compoundList.name(a),srm.compound);
    %ks = find(sn==1);
    ks = strmatch(compoundList.name(a),srm.compound,'exact');
    compoundList.indexMain(a,1) = ks;
    
    h = strmatch(compoundList.name(a),special.cases);

    if isempty(h)
        %these will match as is with no need to dig into the special cases
        %3/12/2021 while I am working on this, just make this broader - match
        %anything with the name and confirm ...
        tn = strcat(compoundList.name(a),' confirm');
        sn = strcmp(tn,srm.compound);
        ks = find(sn==1);
        compoundList.indexConfirm(a,1) = ks;
        clear ks tn sn
    elseif ~isempty(h) && ~isempty(special.match{h})
        %setup the name to use in special - go find that
        s2 = strcmp(special.match(h),srm.compound);
        ks2 = find(s2==1);
        compoundList.indexConfirm(a,1) = ks2;
        clear s ks s2 ks2
    elseif ~isempty(h) && isempty(special.match{h})
        %no confirm ion
        compoundList.indexConfirm(a,1) = NaN;
    end
    clear sn ks h
    
end
clear a


% %%%%% import the CSV file from MAVEN into a dataset array
% %%%%% import the CSV file from MAVEN into a dataset array

%now go ahead and use the dataset function to open the file
importedPeakList = dataset('File',CSVfileFromMAVEN,'delimiter',',');

%clean up
clear CSVfileFromMAVEN tempFile fid tline fidOut tline ans

% warning('off', 'stats:dataset:genvalidnames:ModifiedVarnames');

%%change this if the sampleInfoFile is a CSV file:
sampleInfo = dataset('XLSFile',sampleInfoFile);
% sampleInfo = dataset('File',sampleInfoFile,'delimiter',',','HeaderLines',1) %works for CSV files
clear CSVfileFromMAVEN

%have cases where where files in the TSQ list should not be processed
if pruneData
    k = find(sampleInfo.goodData==0);
    sampleInfo(k,:)=[];
    clear k
end

%MilliQ standard curve or the matrix?
switch curveType
    case 'MQ'
        matchCurve = 'Std Bracket';
    case 'matrix'
        matchCurve = 'Std Bracket matrix';
end

%setup empty matrix for the data and a separate matrix for error
%KL 3/15/2021 updated to pull concentrations from file
switch ionMode
    case 'negative'
        kStandard = find((strcmp(sampleInfo.SampleType,matchCurve)==1) & (strcmp(sampleInfo.ionMode,'negative')==1));           
    case 'positive'
        kStandard = find((strcmp(sampleInfo.SampleType,matchCurve)==1) & (strcmp(sampleInfo.ionMode,'positive')==1));   
    case 'both'
        kStandard = find((strcmp(sampleInfo.SampleType,matchCurve)==1) & (strcmp(sampleInfo.ionMode,'both')==1));          
end

setStandardConcentrations = sampleInfo.concentration_ngML(kStandard); %remember - dataset...not table
        
standardNames = sampleInfo.FileName(kStandard);
nStandards = length(kStandard);

kSample = find(strcmp(sampleInfo.SampleType,'Unknown')==1);
sampleNames = sampleInfo.FileName(kSample);
sampleDetails = sampleInfo.sampleDetails(kSample); %will be more detailed than SampleType
nSamples =length(kSample);
clear kStandard kSample

goodData(1:length(compoundList),nSamples) = NaN; 
goodDataError = goodData;

warning('off', 'stats:dataset:subsasgn:DefaultValuesAddedVariable');
% compoundList.r2_line = 0;
% compoundList.slope = 0;
% compoundList.intercept = 0;
% compoundList.SDslope = 0;
% compoundList.SDintercept = 0;
% compoundList.nPoints = 0;


%Liz wants 5 points in the standard curve (a/o 4/24/2014), so set that here
nRequired = 5;

%skip the labeledStandards - just need those peak areas, change to sending
%in list
% labeledStandards = {'d2-biotin neg','d2-biotin pos'};
exportCurves = table(); %adding 3/24/2021 to send this info outside fxn

for a = 1:length(compoundList);       
    %keep the points for the standard curve - need this outside the function
    exportCurves.mtabName(a) = compoundList.name(a);
    %exportCurves.concentration{a} =  setStandardConcentrations;

    if strcmp(compoundList.name(a),'fumaric acid')
        %stop here for troubleshooting as needed
        %fprintf('here')
        %compoundList.name(a)
    end
    
    clear xdata ydata %here to make checking out one compound easier
    k = strmatch(compoundList.name(a),importedPeakList.compoundId,'exact');
           
    if ~isempty(k)
        %can have cases where nothing good was found. If k is NOT empty,
        %found good data

        smallDS = importedPeakList(k,:); clear k
        try
        [c ia ib] =intersect(smallDS.sample,sampleInfo.FileName);
        catch
            %if you are here - probably bc the wrong file format was
            %exported from El-MAVEN.
            keyboard 
        end
        smallDS.sampleType(ia,1) = sampleInfo.SampleType(ib,1);
        clear c ia ib  

        [~, idxDS, idxStandards] = intersect(smallDS.sample,standardNames);

        %what is the average value in the blanks? Need this for two reasons,
        %(1) to get a zero value for the standard curve and (2) to see if the
        %values in the samples are more/less than what is in the blanks
        kb = find(strcmp(smallDS.sampleType,'Blank')==1);
        %for now, using AreaTop, might play around with that later
        meanBlank = mean(smallDS.peakAreaTop(kb)); clear kb

        xdata = setStandardConcentrations;
        ydata(1:length(xdata),1) = NaN;
        quality = ydata;
        %get all possible values from the standard curve
        ydata(idxStandards) = smallDS.peakAreaTop(idxDS);
        quality(idxStandards) = smallDS.quality(idxDS);

        %%uncomment out next line to export unpruned xdata/ydata
        %exportCurves.show{a} = [xdata ydata quality];
        
        %first, check for linearity (e.g., will drop out points that are 
        %lower at 1000 ng/ml cfd to 500 ng/ml)
                
        %for some compounds, allow lower quality scores...new acids in
        %particular are not so great; KL 4/9/2021
        lowQuality = {'4-methyl-2-oxopentanoic acid',...
            '3-methyl-2-oxopentanoic acid',...
            '2-methyl-4-oxopentanoic acid',...
            '3-3-dimethyl-2oxobutanoic acid',...
            'sucrose_341 neg',...
            'trehalose_387 neg'}';
        if isequal(sum(strcmp(compoundList.name(a),lowQuality)),1)
            kBad = find(quality < 0.05);
            dBad = find(smallDS.quality < 0.05);
        else 
            kBad = find(quality < setQuality);
            dBad = find(smallDS.quality < setQuality);
        end         
            
        ydata(kBad) = NaN;
        xdata(kBad) = NaN; 
        clear quality
        
        %from Winn (10/2016): remove low quality peaks from the data        
        smallDS.peakAreaTop(dBad) = NaN;

        %at this point may have nothing good left...
        if ~isempty(xdata) %only move on if I have something good
            % calculate linearity (see MSdataAnalysis log Jan 2018
            %if <5%, linear ; if >5%, not linear drop point from curve
            setL = 5;
            %cheat...calculate all the curves at once
            for aa = 1:length(xdata) %require 3 points
                dataOut(aa) = getErrors(xdata(1:aa),ydata(1:aa)); %errors for the standard curve

                %see how linear we are...calculate for highest point only
                %do in multiple steps...lost track of ()
                den = (dataOut(aa).slope*xdata(aa))+dataOut(aa).intercept;
                fra = ydata(aa)./den;
                deviation(aa,1) = abs(fra-1)*100;
                clear den fra
            end
            clear aa
            k = find(deviation>setL);
            xdata(k) = NaN;
            ydata(k) = NaN;
            quality(k) = NaN;
            clear k setL deviation


            %MilliQ standard curve or the matrix?
            switch curveType
                case 'MQ'
                    %put the blank at the beginning...remember that this will treat all
                    %blanks equally...this may be a bad thing
                    xdata = cat(1,0,xdata); 
                    ydata = cat(1,meanBlank,ydata);                     
                case 'matrix'
                    %do nothing with meanBlank
            end


            %need to export the xdata and ydata that have been pruned
            exportCurves.xdata{a} = xdata;
            exportCurves.ydata{a} = ydata;
            clear kBad idxDS idxStandards      

            %%%before I go ahead and calculate the values for each sample, make
            %%%sure that the confirm ion is present at the same RT as the main ion.
            %%%adding this in following Melissa's suggestion, 5/19/2014
            i = compoundList.indexConfirm(a);
            %if there is no confirm ion, can't do the check...desthiobiotin neg
            %has no confirm. For now, allow that to pass as if it is always OK
            if ~isnan(i)
                k = strmatch(srm.compound(i),importedPeakList.compoundId,'exact');

                %k will be empty if no confirm ions were found  
                if ~isempty(k) 
                    smallDSconfirm = importedPeakList(k,:); clear k
                    smallDS.cfRT(1:size(smallDS,1),1) = NaN;

                    [c ia ib] = intersect(smallDSconfirm.sample,smallDS.sample);
                    %now compare the retention times for the main ion and the confirm
                    smallDS.cfRT(ib,1) = smallDSconfirm.rt(ia);  
                    %from Winn (10/2016), require the confirm to also meet a 
                    % QC check (for now I am using setQuality/10):
                    %Winn required 0.1, but with the UPLC data this cuts most
                    %thymidine confirm ions
                    %also organize the confirm ion QC data
                    smallDS.cfquality(ib,1) = smallDSconfirm.quality(ia);
                    smallDS.cfsig(ib,1) = smallDSconfirm.signalBaseLineRatio(ia);
                    %set the requirements for the quality of the confirm ion
                    %peak and convert the confirm RT to NaN for those that fall
                    %below the threshold:
                    lowqc = smallDS.cfquality < (setQuality/10) | (smallDS.cfsig <= 1);
                    smallDS.cfRT(lowqc) = NaN;              
                    d = abs(minus(smallDS.cfRT,smallDS.rt));

                    %require that the confirm and the main peak have retention times within
                    %10 seconds (remember to convert to minutes!)
                    %reqRT = 10/60; %this was what I used for the ventDOM project
                    reqRT = 12/60; %MCKS suggested 5/22/2014 that this would be better

                    k2 = find(d > reqRT | isnan(d)); %will be a NaN if peak is not found

                    %for the peaks that fail the confirm retention time check, set them
                    %equal to NaN
                    smallDS.peakAreaTop(k2) = NaN;
                    clear k2 d c ia ib
                elseif isempty(k)
                    %no, will deem this as not good data bc no confirm ions were
                    %found
                    smallDS.peakAreaTop(:) = NaN;
                end
            end

            %remember, will also have cases where no data were found for select samples
            %so need to setup the spacers in there are well
            [c ia ib] = intersect(smallDS.sample,sampleNames);
            tData(1:length(sampleNames),1) = 0; %change to 0, 7/2/2021
            tData(ib) = smallDS.peakAreaTop(ia); 
            clear c ia ib

            %set up a place to keep the data that will be exported
            T = array2table(sampleNames);
            T.sampleDetails = sampleDetails;

            switch curveType
                case 'MQ'
                    %add in two more checks on the data...(1) if the peak area 
                    %is less than the meanBlank, change to a zero 7/2/2021
                    k = find(tData<meanBlank);
                    tData(k) =0;
                    clear k meanBlank                   
                case 'matrix'
                    %do nothing with meanBlank
            end
            


            %another option is to not allow values below the smallest 'good'
            %value in the measured peak areas...add this 5/18/2016
            k = find(tData< min(ydata(2:end)));
            tData(k) = 0; %make zero 7/2/2021
            clear k

            %had been using NaNs to change things with poor quality, change
            %them to 0 here before exporting 7/2/2021
            i = isnan(tData);
            tData(i) = 0; clear t
            T.tData = tData;
            exportCurves.tData{a} = T;
            clear T
        end
        %second check to to if I have good data
    end
end %end of calculating the curves and data
clear a


end %end of getStandardCurves as a function



            


