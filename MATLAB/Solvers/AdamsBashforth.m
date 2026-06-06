function [x, y] = AdamsBashforth(f, a,b, iv, n)

%AdamsBashfort AdamsBashfort's method for solving IVP.
%
%   [x, y] = AdamsBashfort(f,a,b, iv)
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

% Heun as a predictor for the first step

k1 = feval(f,x(1),y(1,:))';
k2 = feval(f,x(2),y(1,:)+h*k1)';
y(2,:) = y(1,:)+h/2*(k1+k2);

% Main program

for k = 2:n-1
    y(k+1,:) = y(k,:)+h/2*(3*feval(f,x(k),y(k,:)) - feval(f,x(k-1),y(k-1,:)))';
end

end