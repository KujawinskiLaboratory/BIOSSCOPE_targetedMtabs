function [calcError, calcConc] = useErrors(myErrorData,measuredSample)
%function [calcError, calcConc] = useErrors(Slope,Intercept,SDslope,SDintercept,measuredSample)
%use the errors on the line to get the errors on the samples actually
%measured and their concentration
%KL 4/21/2014
Intercept = myErrorData.intercept;
Slope = myErrorData.slope;
SDintercept = myErrorData.SDintercept;
SDslope = myErrorData.SDslope;

%calculated concentrations from my hypothetical list
calcConc = (measuredSample - Intercept)./Slope;

%apply to my list of hypothetical unknowns, split this up to make
%it easier to keep track of where the parentheses etc. are
fSQ = (SDintercept./(measuredSample - Intercept)).^2 + (SDslope./Slope).^2;
calcError = calcConc .* sqrt(fSQ);
errorPercent = calcError./calcConc*100;

end %end of useErrors as a function