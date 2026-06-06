function [x, y] = RungeKutta(f, a,b, iv, n)

%RungeKutta Runge-Kutta's method (RK4) for solving IVP.
%
%   [x, y] = RungeKutta(f,a,b, iv)
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
    k2 = feval(f,x(k)+h/2,y(k,:)+h/2*k1)';
    k3 = feval(f,x(k)+h/2,y(k,:)+h/2*k2)';
    k4 = feval(f,x(k+1),y(k,:)+h*k3)';
    y(k+1,:) = y(k,:)+h/6*(k1+2*k2+2*k3+k4);
end

end