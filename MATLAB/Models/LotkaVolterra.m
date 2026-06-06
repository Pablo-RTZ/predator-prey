function f = LotkaVolterra(a,b,c,d)
% LotkaVolterra Returns the ODE for the Lokta-Volterra model
% 
%   f = LotkaVolterra(a,b,c,d)
%
%   Inputs:
%       a - Prey reproduction rate
%       b - Depredation rate
%       c - Predator reproduction rate
%       d - Predator death rate

arguments
a (1,1) double
b (1,1) double
c (1,1) double
d (1,1) double
end


f = @(t, X) [a*X(1)-b*X(1).*X(2); c*X(1).*X(2) - d*X(2)];

end