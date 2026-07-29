function [x, y] = StepsMethod(f, a, b, method, tau, phi, n)

%StepsMethod Solve a DDE using the method of steps.
%
%   [x, y] = StepsMethod(f,a,b,method,tau,phi,n)
%
%   Inputs:
%       f       - DDE right hand side f(t,y,y_delay)
%       a       - Interval start
%       b       - Interval end
%       method  - ODE solver of form method(f,a,b,n)
%       tau     - Delay
%       phi     - History function phi(t)
%       n       - Number of interpolation nodes

arguments
    f (:,1) function_handle
    a (1,1) double
    b (1,1) double {mustBeGreaterThan(b,a)}
    method (1,1) function_handle
    tau (1,1) double {mustBePositive}
    phi (:,1) function_handle
    n (1,1) double {mustBeInteger,mustBePositive} = 10
end

% Initialization

int_n = ceil((b-a) / tau);
n_per = ceil(n / int_n);

% Expanded domain
b_new = a + int_n*tau;

x = [];
y = [];

% Store previous solution for interpolation
x_old = [];
y_old = [];

% Main program

for i = 0:int_n-1

    % Current interval
    a_i = a + i*tau;
    b_i = a + (i+1)*tau;

    % Construct IVP right hand side for this step
    if i == 0
        % Delay comes from history
        f_step = @(t,z) f(t,z,phi(t-tau));
    else
        % Delay uses previously computed solution
        f_step = @(t,z) f(t,z,delay_value(t));
    end

    if i == 0
        iv = phi(a);
    else
        iv = y_old(end,:)';
    end

    % Solve the interval
    [x_i,y_i] = method(f_step,a_i,b_i,iv,n_per);
    x_i = x_i(:);

    if i > 0
        x_i = x_i(2:end);
        y_i = y_i(2:end,:);
    end

    x = [x; x_i];
    y = [y; y_i];

    % Update stored solution
    x_old = x;
    y_old = y;

end


% Nested interpolation function for the delay
    function yd = delay_value(t)

        t_delay = t - tau;

        % Use history if needed
        if t_delay <= a
            yd = phi(t_delay);
        else
            yd = interp1(x_old,y_old,t_delay,'cubic')';
        end

    end

end
