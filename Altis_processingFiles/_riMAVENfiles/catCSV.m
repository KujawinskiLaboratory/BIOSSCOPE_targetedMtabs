function catCSV(folderWithFiles,fileList,useName)
%function catCSV(folderWithFiles,fileList,useName)
%should have done this a while back, make a function to concatenate
%multiple CSV files resulting from MAVEN analysis
%%folderWithFiles --> path (without trailing slash) to files
%%fileList --> the list of files that should be concatenated
%%output will be tempFile which can then be read into MATLAB without issue
%KL 6/25/2021
    tempFile = useName;
    %check if it exists, and if it does throw up an error and stop
    if exist(tempFile)==2;
        fprintf('tempFile exists, deleting it\n')
        delete(tempFile)
    end

    %open up the first CSV file   
    fid1 = fopen([folderWithFiles filesep fileList{1}]);

    tline = fgetl(fid1);
    if ~strcmp(tline(end),',')
        tline = strcat(tline,','); %put the comma at the end of the headerline
    end
    %open up the temp file, and set it to be write-able
    fidOut = fopen(tempFile,'w');
    fprintf(fidOut, '%s\n',tline); %write the first line

    tline = fgetl(fid1); %get the next line of the file
    while ischar(tline)
        fprintf(fidOut, '%s\n', tline);
        tline = fgetl(fid1);    
    end
    fclose(fid1); %close the first file
    clear fid1 tline
    
    if length(fileList) > 1
        % pull in the rest of the files using a loop
        for ai = 2:length(fileList)
            %fprintf('In loop %d of %d\n',a,length({fileList}))
            %now get the next file...since the headerline will be the same can skip it
            fid2 = fopen([folderWithFiles filesep fileList{ai}]);
            tline = fgetl(fid2);
            tline = fgetl(fid2); %get the next line of the file
            while ischar(tline)
                fprintf(fidOut, '%s\n', tline);
                tline = fgetl(fid2);    
            end
            fclose(fid2);
            clear fid2 tline
        end
        clear ai
    end

    fclose(fidOut);
    clear D ans fidOut wDir
end
