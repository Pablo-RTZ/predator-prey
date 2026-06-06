function f = FunctionalResponse(a,b,c,d,h,K)
% FunctionalRespone Models the predator-prey problem with functional
% response (handling time)
% 
%   f = FunctionalResponse(a,b,c,d)
%
%   Inputs:
%       a - Prey reproduction rate
%       b - Depredation rate
%       c - Predator reproduction rate
%       d - Predator death rate
%       h - Handling time
%       K - Prey carrying capacity

arguments
a (1,1) double
b (1,1) double
c (1,1) double
d (1,1) double
h (1,1) double
K (1,1) double
end

f = @(t, X) [a*X(1).*(1-X(1)/K)-b*X(1).*X(2)./(1+b*h*X(1)); c*b*X(1).*X(2)./(1+b*h*X(1)) - d*X(2)];

end