function f = DelayedLotkaVolterra(a,b,c,d)
%DelayedLogLotkaVolterra Returns the DDE for the Lotka-Volterra model
% with delayed predator gestation.
%
%   f = DelayedLogLotkaVolterra(a,b,c,d,)
%
%   Inputs:
%       a - Prey reproduction rate
%       b - Depredation rate
%       c - Predator reproduction rate
%       d - Predator death rate
%
%   Output:
%       f - Function handle f(t,X,X_tau)

arguments
    a (1,1) double
    b (1,1) double
    c (1,1) double
    d (1,1) double
end

f = @(t,X,X_tau) [
    a*X(1) - b*X(1).*X(2);
    c*X_tau(1).*X_tau(2) - d*X(2)
];

end
