function [t, P,Z] = EulerDiffusion(f, g, P0, Z0, Dp, Dz, t_end, n)
%EulerDiffusion Euler's method for solving reaction-diffusion problems.
%
%   [t, P, Z] = EulerDiffusion(f, g, P0, Z0, Dp, Dz, t_end, n)
%   Solves the problem and returns function at evaluation nodes
%
%   Inputs:
%       f       - Function handle for the prey reaction
%       g       - Function handle for the predator reaction
%       P0      - Initial condition for P
%       Z0      - Initial condition for Z
%       Dp      - Diffusion coeficient for P
%       Z0      - Diffusion coeficient for Z
%       t_end   - Time interval end
%       n       - Number of timesteps

arguments
    f (1,1) function_handle
    g (1,1) function_handle
    P0 (:,:) double
    Z0 (:,:) double
    Dp (1,1) double
    Dz (1,1) double
    t_end (1,1) double
    n (1,1) double {mustBeInteger,mustBePositive}

end

if ~isequal(size(P0), size(Z0))
    error('Initial conditions must have the same size.');
end

P(:,:,1) = P0;
Z(:,:,1) = Z0;
[Nx, Ny] = size(P0);

dx = 1/(Nx-1); dy = 1/(Ny-1);
x = linspace(0,1,Nx); y = linspace(0,1,Ny);

t = linspace(0,t_end,n);
dt = t_end/(n-1);

for i = 1:n-1

    Pk = P(:,:,i);
    Zk = Z(:,:,i);

    % Laplacian of P
    LP = zeros(Nx,Ny);

    % interior
    LP(2:end-1,2:end-1) = ...
        (Pk(3:end,2:end-1) - 2*Pk(2:end-1,2:end-1) + Pk(1:end-2,2:end-1))/dx^2 + ...
        (Pk(2:end-1,3:end) - 2*Pk(2:end-1,2:end-1) + Pk(2:end-1,1:end-2))/dy^2;

    % Neumann boundaries
    LP(1,:)   = 2*(Pk(2,:)     - Pk(1,:))/dx^2;
    LP(end,:) = 2*(Pk(end-1,:) - Pk(end,:))/dx^2;
    LP(:,1)   = LP(:,1)   + 2*(Pk(:,2)     - Pk(:,1))/dy^2;
    LP(:,end) = LP(:,end) + 2*(Pk(:,end-1) - Pk(:,end))/dy^2;

    % Laplacian of Z
    LZ = zeros(Nx,Ny);

    % interior
    LZ(2:end-1,2:end-1) = ...
        (Zk(3:end,2:end-1) - 2*Zk(2:end-1,2:end-1) + Zk(1:end-2,2:end-1))/dx^2 + ...
        (Zk(2:end-1,3:end) - 2*Zk(2:end-1,2:end-1) + Zk(2:end-1,1:end-2))/dy^2;

    % Neumann boundaries
    LZ(1,:)   = 2*(Zk(2,:)     - Zk(1,:))/dx^2;
    LZ(end,:) = 2*(Zk(end-1,:) - Zk(end,:))/dx^2;
    LZ(:,1)   = LZ(:,1)   + 2*(Zk(:,2)     - Zk(:,1))/dy^2;
    LZ(:,end) = LZ(:,end) + 2*(Zk(:,end-1) - Zk(:,end))/dy^2;

    % Explicit Euler update
    P(:,:,i+1) = P(:,:,i) + dt*((f(t(i),P(:,:,i),Z(:,:,i)))+ Dp*LP);
    Z(:,:,i+1) = Z(:,:,i) + dt*((g(t(i),P(:,:,i),Z(:,:,i)))+ Dz*LZ);
end

end