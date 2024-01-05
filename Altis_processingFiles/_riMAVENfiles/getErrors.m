function dataOut = getErrors(xdata,ydata);
%function dataOut = getErrors(xdata,ydata);
%From this web site:
%http://terpconnect.umd.edu/~toh/spectrum/LeastSquaresMatlab.txt
%KL modifying 4/21/2014
%I don't think I need this in here, but leave for now KL 7/1/2021

x = xdata;
y = ydata;
% Simple Matlab script for calculating the first-order least-square fit of y vs x,
% including the Slope and Intercept and the predicted standard deviation of
% the slope (SDSlope) and intercept (SDIntercept).

NumPoints=length(x);
Sxx = sum((x-mean(x)).^2);
Syy = sum((y-mean(y)).^2);
Sxy = sum((x-mean(x)).*(y-mean(y)));
Slope = Sxy./Sxx;
Intercept = mean(y)-Slope*mean(x);

Sy = sqrt((Syy-Slope^2*Sxx)/(NumPoints-2));

SDslope = Sy/sqrt(Sxx);
SDintercept = Sy*sqrt(1./(NumPoints-(sum(x).^2)./sum(x.^2)));

r2 = 1 - ((Syy-Slope^2*Sxx) ./Syy);

%data to send out of this function (when it is a function)
dataOut.slope = Slope;
dataOut.intercept = Intercept;
dataOut.SDslope = SDslope;
dataOut.PercentSlopeError = SDslope./Slope;
dataOut.SDintercept = SDintercept;
dataOut.PercentInterceptError = SDintercept./Intercept;
dataOut.r2 = r2;

end %end of getErrors as a function