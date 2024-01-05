function dataOut = considerMAVEN_stepC(xdata,ydata,maxInSamples)
%function dataOut = considerMAVEN_stepC(xdata,ydata,maxInSamples)
%need xdata & ydata for the curve itself 
%maxInSamples --> what is the highest point on the standard curve?
%send out the slope/intercept etc. from the standard curve
%KL 7/2/2021

%Liz wants 5 points in the standard curve (a/o 4/24/2014), so set that here
nRequired = 5;

%need to deal with the idea of how big to allow the curve to be
%and, what is the max value in my samples? should have at least one point above that
%if all the unknowns fail the quality check, this next step will fail. 
if isnan(maxInSamples)
    %easiest to make kMax empty
    kMax = [];
else
    kMax = find(ydata <= maxInSamples);
end

if ~isempty(kMax)
    %have at least one point on the curve
    if isequal(kMax(end),length(ydata)) 
        %already at the end of the standard curve...so use all the points
        %do nothing...but send up a flag since the data are above
        %the standard curve
        %disp([compoundList.name(a) ' is above the standard curve'])
        %fprintf('here')
    elseif isequal(kMax(end)+1,length(ydata));
        %only one more above the points in the standard curve, use all the
        %points
    elseif isequal(kMax,1);
        %data are at the low end of the standard curve, but let's require 
        %more points above my data to get a reasonable curve...
        xdata = xdata(1:nRequired);
        ydata = ydata(1:nRequired);
    elseif length(kMax)+2  < nRequired
        % use the number of points sent in nRequired
        ydata = ydata(1:nRequired);
        xdata = xdata(1:nRequired);
    elseif length(kMax) + 1 < nRequired
        %use the standard curve to one point beyond the range of my
        %samples
        ydata = ydata(1:kMax(end)+1);
        xdata = xdata(1:kMax(end)+1);
    else
        %use the standard curve to one point beyond the range of my
        %samples
        ydata = ydata(1:kMax(end)+1);
        xdata = xdata(1:kMax(end)+1);

    end
elseif isempty(kMax)
    %all of the points in the standard curve are higher than what was
    %measured in the samples; can truncate to the number of points required
    ydata = ydata(1:nRequired);
    xdata = xdata(1:nRequired);

end
clear kMax maxInSamples

%need at least three points to make a curve AND get the error estimates
try
    show = [xdata ydata];
catch 
    error('Something wrong, xdata and ydata should be the same size')
end

i = isnan(show);
sfmi = sum(i,2);
k = find(sfmi==0);
xdata = xdata(k);
ydata = ydata(k); 
clear show i sfmi k

if length(xdata)>2 
        %I have enough for a standard curve; getErrors is now its own
        %function in the _riMAVENfiles folder (keep the name, getErrors was
        %within considerMAVEN in its previous rendition)
        dataOut = getErrors(xdata,ydata); %errors for the standard curve
        dataOut.nPoints = length(ydata); 
elseif length(xdata)<=2
        %Get here when there is no curve passing the quality checks
        dataOut = {'QCfail'};
end

end %end of the calculateCurve function

    