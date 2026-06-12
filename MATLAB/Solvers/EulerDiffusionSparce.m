function [t, P,Z] = EulerDiffusionSparce(f, g, P0, Z0, Dp, Dz, t_end, n)
%EulerDiffusionSparce Euler's method for solving reaction-diffusion
% problems using sparce matrices.
%
%   [t, P, Z] = EulerDiffusionSparce(f, g, P0, Z0, Dp, Dz, t_end, n)
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

[Nx, Ny] = size(P0);

dx = 1/(Nx-1); dy = 1/(Ny-1);

t = linspace(0,t_end,n);
dt = t_end/(n-1);

% 1D Neumann Laplacian in x
ex = ones(Nx,1);
Lx = spdiags([ex -2*ex ex],[-1 0 1],Nx,Nx);

% Neumann BCs
Lx(1,1)     = -2;
Lx(1,2)     =  2;
Lx(end,end) = -2;
Lx(end,end-1)= 2;

Lx = Lx/dx^2;

% 1D Neumann Laplacian in y
ey = ones(Ny,1);
Ly = spdiags([ey -2*ey ey],[-1 0 1],Ny,Ny);

% Neumann BCs
Ly(1,1)      = -2;
Ly(1,2)      =  2;
Ly(end,end)  = -2;
Ly(end,end-1)=  2;

Ly = Ly/dy^2;

% 2D Laplacian
Ix = speye(Nx);
Iy = speye(Ny);

L = kron(Iy,Lx) + kron(Ly,Ix);   % sparse (Nx*Ny)-by-(Nx*Ny)

P = zeros(Nx*Ny,n);
Z = zeros(Nx*Ny,n);

P(:,1) = reshape(P0,[],1);
Z(:,1) = reshape(Z0,[],1);

for k = 1:n-1

    Pk = reshape(P(:,k),Nx,Ny);
    Zk = reshape(Z(:,k),Nx,Ny);

    reactionP = f(t(k),Pk,Zk);
    reactionZ = g(t(k),Pk,Zk);

    P(:,k+1) = P(:,k) + dt*( ...
        reshape(reactionP,[],1) + Dp*(L*P(:,k)));

    Z(:,k+1) = Z(:,k) + dt*( ...
        reshape(reactionZ,[],1) + Dz*(L*Z(:,k)));

end

P = reshape(P, Nx, Ny, n);
Z = reshape(Z, Nx, Ny, n);

end