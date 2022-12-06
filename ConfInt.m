%Confidence Interval for Linear regression:

function [yconf] = ConfInt(x,y,yhat,p,type)
%input: x and y datapoints, yhat = fit line evaluated at x, p = confidence
%(.05 for 95%)

%degrees of freedom:
n = length(x);
df = n-2;
% Calculate residuals:
res = y - yhat; ress = res.^2;

%residual standard error:
sy = (sum(ress)/df).^0.5;


%two - sided t distribution for 1-p
tstar = tinv(1-p/2,df);

%standard deviation of x:
sx = std(x);

if strcmp(type,'Confidence')
    %calculate upper y bounds:
    yconf(:,1) = yhat + tstar*sy*(1/n + (x-mean(x)).^2/((n-1)*sx^2)).^0.5;
    yconf(:,2) = yhat - tstar*sy*(1/n + (x-mean(x)).^2/((n-1)*sx^2)).^0.5;
elseif strcmp(type,'Prediction')
    %calculate upper y bounds:
    yconf(:,1) = yhat + tstar*sy*(1+ 1/n + (x-mean(x)).^2/((n-1)*sx^2)).^0.5;
    yconf(:,2) = yhat - tstar*sy*(1+ 1/n + (x-mean(x)).^2/((n-1)*sx^2)).^0.5;  
end
%put in order for plotting:
% yconf = sort(yconf);
% xcord = sort(x);
end

