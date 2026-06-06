function [x, y] = DormandPrince(f, a,b, iv, n)

%DormandPrince Dormand-Prince method (DOPRI5) for solving IVP.
%
%   [x, y] = DormandPrince(f,a,b, iv)
%   Solves the problem and returns function at evaluation nodes
%
%   Inputs:
%       f   - ODE
%       a   - Interval start
%       b   - Interval end
%       iv  - Initial value
%       n   - Number of nodes to interpolate

arguments
    f (:,1) function_handle
    a (1,1) double
    b (1,1) double {mustBeGreaterThan(b,a)}
    iv (:,1) double
    n (1,1) double {mustBeInteger,mustBePositive} = 10
end

% Initialization

h = (b-a)/n;
x = a:h:b;
y = zeros(n+1,length(iv));
y(1,:) = iv;

% Main program

for k = 1:n

    k1 = feval(f,x(k),y(k,:))';

    k2 = feval(f,...
        x(k) + h*(1/5),...
        y(k,:) + h*(1/5)*k1)';

    k3 = feval(f,...
        x(k) + h*(3/10),...
        y(k,:) + h*(3/40*k1 + 9/40*k2))';

    k4 = feval(f,...
        x(k) + h*(4/5),...
        y(k,:) + h*(44/45*k1 - 56/15*k2 + 32/9*k3))';

    k5 = feval(f,...
        x(k) + h*(8/9),...
        y(k,:) + h*(19372/6561*k1 - 25360/2187*k2 + ...
                    64448/6561*k3 - 212/729*k4))';

    k6 = feval(f,...
        x(k) + h,...
        y(k,:) + h*(9017/3168*k1 - 355/33*k2 + ...
                    46732/5247*k3 + 49/176*k4 - ...
                    5103/18656*k5))';

    k7 = feval(f,...
        x(k) + h,...
        y(k,:) + h*(35/384*k1 + 500/1113*k3 + ...
                    125/192*k4 - 2187/6784*k5 + ...
                    11/84*k6))';

    % 5th-order Dormand-Prince solution
    y(k+1,:) = y(k,:) + h*( ...
          35/384*k1 ...
        + 500/1113*k3 ...
        + 125/192*k4 ...
        - 2187/6784*k5 ...
        + 11/84*k6 );

end

end