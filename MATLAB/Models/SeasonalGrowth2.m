function f = SeasonalGrowth2(A,b,c,d,h,K,alpha,phi)
% SeasonalGrowth2 Models the predator-prey problem with functional
% response and seasonal growth and death rate
% 
%   f = SeasonalGrowth2(A,b,c,d,h,K,alpha,phi)
%
%   Inputs:
%       A     - Average prey reproduction rate
%       b     - Depredation rate
%       c     - Predator reproduction rate
%       d     - Predator death rate
%       h     - Handling time
%       K     - Prey carrying capacity
%       alpha - Predator death rate amplitude
%       phi   - Predator death rate phase (months)

arguments
A (1,1) double
b (1,1) double
c (1,1) double
d (1,1) double
h (1,1) double
K (1,1) double
alpha (1,1) double
phi (1,1) double
end

f = @(t, X) [A*(1/2+cos(pi*t/12).^2).*X(1).*(1-X(1)/K)-b*X(1).*X(2)./(1+b*h*X(1));
    c*b*X(1).*X(2)./(1+b*h*X(1)) - d*(1+alpha*cos(pi*(t+phi)/6)).*X(2)];

end