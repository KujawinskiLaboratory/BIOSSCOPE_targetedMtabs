function newOne = stripName(one)
%function newOne = stripName(one)
%the metabolite names are really tedious, so I am going to strip all the
%annoying pieces out so that the matching is easier. Do this as a function
%and then set it loose everyone
%KL 10/5/2022
%There are three steps needed to make these names manageable
%1. strip out all the apostrophes --> char(39)
%1b. strip out the hyphens --> char(45)
%2. make names all lower
%3. deblank the names
if ischar(one)  
    %char will only be one metabolite
    %do nothing...already ready for next step
    one(one==char(39))=[];
    one(one==char(45))=[];
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
