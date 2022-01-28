function [sampleNames, keepGoodData] = considerMAVEN_stepB_v2(tableOne,tableMix8,labeledCpds)
%function [sampleNames, keepGoodData] = considerMAVEN_stepB_v2(tableOne,tableMix8,labeledCpds)
%tableOne --> table with information about all the compounds in one ion
%mode (either positive or negative)
%tableMix8 --> table with the information about the mix8 standard curves
%which have been run separately and include a combination of pos and neg
%mode metabolites
%labeledCpds --> the list of isotopically-labeled compounds to be 
%given special treatment (no concentrations, just keep peak areas)
%parts of this were originally in considerMAVEN
%KL 6/28/2021 and after all the work to add a third set for Mix 8, we go
%back to two SRM lists, but need three standard curves (pos/neg/Mix8). This
%is run after getStandardCurves_v1.m 

% warning('off','MATLAB:table:ModifiedAndSavedVarnames')
% warning('off','stats:dataset:ModifiedVarnames')
%warning('off')

%go through each compound listed in tableOne
   
%need to figure out what the range of samples is BEFORE calculating the
%standard curve. So, pull those pieces into here.
%set up a column to mark stuff
tableOne.flag(:) ={[]}; %dumb MATLAB option for adding empty column

for a = 1:size(tableOne,1);
%     %stop here for troubleshooting as needed   
%     if strcmp(tableOne.mtabName(a),'trehalose_387 neg')
%         %fprintf('here')
%         %tableOne.mtabName(a)
%     end
    
    %sampleNames will have been repeated in tableOne (oops) check the order
    if ~isempty(tableOne.tData{a})
        sampleNames = tableOne.tData{a}.sampleNames;
    end
        
    if ~isempty(tableOne.tData{a}) & ~isequal(tableOne.tData{a}.sampleNames,sampleNames)
        error('something is wrong, these should match \n')
    end
    
    %four possibilities...different path for each    
    %set up to look for (1) labeledCpd, (2) mix8 curve, (3) pos or neg
    %curve, (4) or no data at all
    if isempty(tableOne.xdata{a}) & ~contains(tableOne.mtabName(a),labeledCpds)
        %do nothing, no data are available (option #4 on list)
    elseif sum(strcmp(tableOne.mtabName(a),labeledCpds))>0 %option #1 on list
        %this is a labeled compound, just keep the peak areas, (IF they exist)
        if ~isempty(tableOne.xdata{a})
            tableOne.goodData(a,:) = tableOne.tData{a}.tData';
        end
    else %have data to use in some manner (option #2 or #3)
        if sum(strcmp(tableOne.mtabName(a),tableMix8.mtabName))>0 %option #2: need mix8 for curve
            %use the standard curve from Mix8, first have to go find it.
            s = strcmp(tableOne.mtabName{a},tableMix8.mtabName);
            ks = find(s==1);
            
            %use that index to get the info on standard curve from the mix 8 run
            try
            xdata = tableMix8.xdata{ks};
            catch
                fprintf('here')
            end
            
            ydata = tableMix8.ydata{ks};
            clear s ks
            
            %tData (peaks areas of samples) still comes from tableOne
            tData = tableOne.tData{a}.tData;
            
            %figure out which are the unknowns in tData
            s = strcmp(tableOne.tData{a}.sampleDetails,'Unknown');
            ks = find(s==1);
            %and then find the max value in samples
            maxInSamples = max(tData(ks)); 
            clear s ks
            
            %calculate parameters of the standard curve, if length(xdata>2)
            dataOut = considerMAVEN_stepC_v1(xdata,ydata,maxInSamples);

            if strcmp(dataOut,'QCfail')
                 calcConc = [];
                 %put in a marker
                  tableOne.flag(a) = {'QCfail'};

            else
                %then, calculate the concentrations in samples
                %useErrors is function that was originally within considerMAVEN
                %keep this name 
                [calcError, calcConc] = useErrors(dataOut,tData);   
            end
        else % this curve is here in pos/neg data; use it
            %use that index to get the info on standard curve from the mix 8 run
            xdata = tableOne.xdata{a};
            ydata = tableOne.ydata{a};
            tData = tableOne.tData{a}.tData;
            
            %figure out which are the unknowns in tData
            s = strcmp(tableOne.tData{a}.sampleDetails,'Unknown');
            ks = find(s==1);
            %and then find the max value in Unknowns
            maxInSamples = max(tData(ks)); 
            clear s ks
            
            dataOut = considerMAVEN_stepC_v1(xdata,ydata,maxInSamples);
            if strcmp(dataOut,'QCfail');
                calcConc = [];
                tableOne.flag(a) = {'QCfail'};
            else
                [calcError, calcConc] = useErrors(dataOut,tData);
            end
        end
        
        %now that I have the data, clean it up before stowing it away, just
        %put away good data, otherwise will be all zeros
        if ~isempty(calcConc)
            %add some check, first, if the slope is negative, the whole set is garbage
            if dataOut.slope <0 
                calcConc(:) = 0;
                calcError(:) = 0;
            end

            kz = find(calcConc < 0); %easiest: sample is less than zero
            calcConc(kz) = 0;
            clear kz

            kz = find(calcConc - calcError < 0);
            calcConc(kz) = 0;
            clear kz

            %now stuff the data back into the table
            tableOne.goodData(a,:) = calcConc';
            tableOne.goodDataError(a,:) = calcError';
        end
        
        clear calcConc calcError
    end
            
       
end % end of the list of compounds to go through
clear a

%keep a variant of this syntax to make downstream life the same. Remember
%that kgd is a dataset...change that later.
keepGoodData = table2dataset(tableOne);
keepGoodData.name = keepGoodData.mtabName; %duplicate for now

%leave in zeros, but delete metabolites where the curve fails the quality
%check as I could not measure those if I wanted to. These will either be
%marked with QCfail or will be empty
s = strcmp(keepGoodData.flag,'QCfail');
ks = find(s==1);
keepGoodData(ks,:)=[];
clear s ks


end %end of calculateConcentrations_v2.m as a function



            


