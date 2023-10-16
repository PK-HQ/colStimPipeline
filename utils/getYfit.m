function [yFitValue,coefficients]=getYfit(x,y,fitExponent,desired_xvalue)
coefficients = polyfit(x,y,fitExponent);
yFitValue = coefficients(1).*desired_xvalue + coefficients(2);
end