function f = Logistic(a,b,c,d, K)
% Logistic Returns the ODE for the Lokta-Volterra model with logistic
% growth
% 
%   f = Logistic(a,b,c,d)
%
%   Inputs:
%       a - Prey reproduction rate
%       b - Depredation rate
%       c - Predator reproduction rate
%       d - Predator death rate
%       K - Carrying capacity

arguments
a (1,1) double
b (1,1) double
c (1,1) double
d (1,1) double
K (1,1) double
end


f = @(t, X) [a*X(1).*(1- X(1)/K)-b*X(1).*X(2); c*X(1).*X(2) - d*X(2)];

end