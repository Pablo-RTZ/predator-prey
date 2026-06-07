function f = SeasonalGrowth(A,b,c,d,h,K)
% SeasonalGrowth Models the predator-prey problem with functional
% response and seasonality
% 
%   f = SeasonalGrowth(A,b,c,d)
%
%   Inputs:
%       A - Average prey reproduction rate
%       b - Depredation rate
%       c - Predator reproduction rate
%       d - Predator death rate
%       h - Handling time
%       K - Prey carrying capacity

arguments
A (1,1) double
b (1,1) double
c (1,1) double
d (1,1) double
h (1,1) double
K (1,1) double
end

f = @(t, X) [A*(1/2+cos(pi*t/12).^2).*X(1).*(1-X(1)/K)-b*X(1).*X(2)./(1+b*h*X(1)); c*b*X(1).*X(2)./(1+b*h*X(1)) - d*X(2)];

end