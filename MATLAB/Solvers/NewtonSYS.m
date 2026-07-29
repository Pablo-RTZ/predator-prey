function [sol,dif,res,x,iter, ACOC] = NewtonSYS(F, dF, x0, opts)

%NewtonSYS Newton-Raphson's method for solving nonlinear systems.
% Uses function handles in vector form.
%
%   sol = NewtonSYS(@F, @dF, x0)
%   Uses a vector of 0s as initial guess, a 1e-10 tolerance, and 50 iterations by
%   default
%
%   sol = NewtonSYS(@F, @dF, x0, "tol", 1e-8, "maxiter", 100)
%   allows name-value pair inputs in any order.
%
%   Outputs:
%       sol   - root approximation
%       dif   - step sizes (vector)
%       res - function evaluation at nodes |F(x0)|
%       x     - solution aproximations
%       iter  - number of iterations
%       ACOC  - Approximate Computational Order of Convergence

arguments
    F (:,1) function_handle
    dF (:,:) function_handle
    x0 (:,1) double
    opts.tol (1,1) double = 1e-10
    opts.maxiter (1,1) double {mustBeInteger, mustBeNonnegative} = 50
end

% Initialization

iter = 0;
dif = opts.tol+1;
res = opts.tol+1;

F_x0 = F(x0);
dF_x0 = dF(x0);

% Main program

while iter < opts.maxiter && dif(end)+res(end) > opts.tol
    iter = iter+1;

    x(:,iter) = x0 - (dF_x0 \ F_x0);

    dif(iter) = norm(x(:,iter) - x0);

    x0 = x(:,iter);
    F_x0 = F(x0);
    dF_x0 = dF(x0);
    res(iter) = norm(F_x0);

end

if iter < opts.maxiter
    sol = x(:,end);
    ACOC = log(res(3:end)./res(2:end-1))./log(res(2:end-1)./res(1:end-2));
else
    sol = NaN;
    ACOC = NaN;
end

% Optional arguments clearing

if nargout < 6, clear ACOC; end
if nargout < 5, clear iter; end
if nargout < 4, clear x; end
if nargout < 3, clear res; end
if nargout < 2, clear dif; end

end